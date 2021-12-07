# Title     : Limma analysis
# Objective : Mostly to take advantage of the duplicateCorrelation() function in Limma
# Created by: marcus
# Created on: 07/02/2020

LimmaTesting <- R6::R6Class("LimmaTesting", list(
	pathing = NULL,
	sample_details = NULL,
	sample_details_source = NULL,
  	reference_details = NULL,
	t2g = NULL,
	dge_list = NULL,
	txi_kallisto = NULL,
	id = NULL,
	initialize = function(pathing, sample_details, reference_details, id = "sample") {
		if (!id %in% colnames(sample_details$s2c)) {
			print(paste(id, "is not found in the sample details file. Consider supplying a custom value for id. Aborting"))
			stop()
		}
	  	self$pathing <- pathing
	  	self$sample_details_source <- sample_details
		self$sample_details <- sample_details$clone()
		self$reference_details <- reference_details
	  	self$id <- id
		return(invisible(self))
	},
	tx_import = function(txi_kallisto = NULL) {
	  	if (!(is.null(txi_kallisto))) {
		 	self$txi_kallisto <- txi_kallisto
		} else {
			files <- file.path(self$sample_details_source$s2c$path, "abundance.h5")
			self$txi_kallisto <- tximport::tximport(
				files,
				type = 'kallisto',
				tx2gene = self$reference_details$t2g,
				# From the bioconductor package notes,
				# "Because limma-voom does not use the offset matrix stored in y$offset,
				#   we recommend using the scaled counts generated from abundances,
				#   either "scaledTPM" or "lengthScaledTPM"
				countsFromAbundance = "lengthScaledTPM"
			)
		}
		return(invisible(self))
	},
	prepare_dgelist = function(group = NULL) {
	  	if (!(group %in% colnames(self$sample_details$s2c))) {
	  		stop('Grouping variable not found in sample details')
	  	}
	 	self$dge_list <- edgeR::DGEList(
		  counts = self$txi_kallisto$counts,
		  # here because we are looking at genes and not transcript isoforms
			# I just take the mean estimate of effective length rather than considering any length bias
		  genes = data.frame(rowMeans(self$txi_kallisto$length)),
		  samples = self$sample_details_source$s2c,
		  group = self$sample_details_source$s2c[, group]
		)
		if (length(self$sample_details_source$s2c[,self$id]) != length(unique(self$sample_details_source$s2c[,self$id]))) {
			print("Technical replicates will be summed.")
			self$dge_list <- edgeR::sumTechReps(self$dge_list, self$sample_details_source$s2c[, self$id])
			all_same <- function(v) {
			  if(length(unique(v)) == 1) {
				TRUE
			  } else {
				FALSE
			  }
			}
			s2c_same <- self$sample_details_source$s2c %>%
			  dplyr::group_by(across(all_of(self$id))) %>%
			  dplyr::summarise(across(.fns=all_same))
			keep_columns <- unlist(apply(s2c_same, 2, all))
			keep_columns[self$id] <- T
			self$sample_details$s2c <- self$sample_details_source$s2c %>%
			  dplyr::group_by(across(all_of(self$id))) %>%
			  dplyr::slice(1) %>%
			  dplyr::select(which(keep_columns)) %>%
			  dplyr::ungroup()
		}
	  	# self$dge_list <- edgeR::calcNormFactors(self$dge_list)
	  	# best to not filter on normalised counts as TMM normalisation fails at low counts
	    # see comment by Aaron Lun: https://support.bioconductor.org/p/116351/
	   	self$dge_list <- self$dge_list[edgeR::filterByExpr(self$dge_list), ] %>%
	  		edgeR::calcNormFactors() # renormalise after filter
	    normalised_counts <- as.data.frame(edgeR::cpm(self$dge_list$counts)) %>%
		  tibble::rownames_to_column("GeneID")
		write.table(
		  normalised_counts,
		  file=file.path(pathing$base_path, paste0("cpm_filtered_across",fs::path_sanitize(self$id),".txt")),
		  sep="\t",
		  row.names = F
		)
	},
	design = function(model) {
	  	design <- model.matrix(model, self$sample_details$s2c)
	  	return(design)
	},
	calculate_correlation =function (voom_object, design, block) {
		corfit <- limma::duplicateCorrelation(
		  voom_object,
		  design,
		  block = block
		)
	  	print(paste('Estimated average correlation between block members:', corfit$cor))
	  	return(corfit$cor)
	},
	fit = function(model, block = NULL) {
	  	design <- self$design(model)
	  	voom_object <- limma::voom(self$dge_list, design)
	  	if (is.null(block)) {
		  	fit <- limma::lmFit(voom_object, design)
		} else {
		  	blocking_vector <- self$sample_details$s2c[, block]
		  	cor <- self$calculate_correlation(voom_object, design, blocking_vector)
		  	fit <- limma::lmFit(
			  voom_object,
			  block = blocking_vector,
			  correlation = cor
			)
		}
	  	fit <- limma::eBayes(fit)
		return(fit)
	},
	contrast_fit = function(model, contrasts, block = NULL) {
	  fit <- self$fit(model, block)
	  contrast_fit <- limma::eBayes(limma::contrasts.fit(fit, contrasts))
	  return(contrast_fit)
	},
	deg_results = function(fit, coef_index, n = NULL) {
	  	# coeff_index is the indexed position of the cofficient of interest in the model,
	  	# i.e. which collumn in model.matrix
	  	if (is.null(n)) {
		  n <- length(fit$F)
		}
	  	limma::topTable(fit, coef = coef_index, number = n) %>%
		   tibble::rownames_to_column('GeneID')
	},
	ndeg = function(fit) {
	  	summary(limma::decideTests(fit))
	},
	plot_volcano = function (
		fit,
		coef_index,
		file_path,
		top_n = 10,
		point_alpha = 0.2
	) {
		res <- self$deg_results(fit, coef_index) %>% dplyr::mutate(
			.,
			significant = as.factor(adj.P.Val <= 0.05),
		)
		colours_vector <- c('black', 'red')
		p <- ggplot2::ggplot(res, ggplot2::aes(logFC, -log10(adj.P.Val), label=GeneID)) +
			ggplot2::geom_point(ggplot2::aes(colour = significant,), alpha = point_alpha) +
			ggrepel::geom_text_repel(
				data = res[1:top_n,],
				min.segment.length = ggplot2::unit(0,'lines')
			) +
			ggplot2::scale_colour_manual(values = colours_vector) +
			ggplot2::xlab("LogFC") +
			ggplot2::ylab("-log10(adj.P.Val)") +
			ggplot2::geom_vline(xintercept = 0, colour = "black", linetype = "longdash")
		ggplot2::ggsave(file_path, plot=p)
	}

))