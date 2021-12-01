# Title     : Limma analysis
# Objective : Mostly to take advantage of the duplicateCorrelation() function in Limma
# Created by: marcus
# Created on: 07/02/2020

LimmaTesting <- R6::R6Class("LimmaTesting", list(
	base_path = NULL,
	sample_details = NULL,
  	reference_details = NULL,
	t2g = NULL,
	dge_list = NULL,
	txi_kallisto = NULL,
	outliers = NULL,
	initialize = function(pathing, sample_details, reference_details, outliers = NULL) {
		self$base_path <- pathing$return_to_base()$append_to_path('quality_control')$current_path
		self$sample_details <- sample_details
		self$reference_details <- reference_details
	  	self$outliers <- outliers
		return(invisible(self))
	},
	tx_import = function(txi_kallisto = NULL) {
	  	if (!(is.null(txi_kallisto))) {
		 	self$txi_kallisto <- txi_kallisto
		} else {
			files <- file.path(self$sample_details$s2c$path, "abundance.h5")
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
		colnames(self$txi_kallisto$counts) <- self$sample_details$s2c$sample
		return(invisible(self))
	},
	prepare_dgelist = function(group, ID=NULL) {
	 	self$dge_list <- edgeR::DGEList(
		  counts = self$txi_kallisto$counts[, !(colnames(self$txi_kallisto$counts) %in% self$outliers)],
		  # here because we are looking at genes and not transcript isoforms
			# I just take the mean estimate of effective length rather than considering any length bias
		  genes = data.frame(rowMeans(self$txi_kallisto$length)),
		  samples = self$sample_details$s2c[!(self$sample_details$s2c$sample %in% self$outliers),],
		  group = self$sample_details$s2c[!(self$sample_details$s2c$sample %in% self$outliers), colnames(self$sample_details$s2c) == group]
		)
	  	if (!(is.null(ID))) {
		  self$dge_list <- edgeR::sumTechReps(lt$dge_list, self$sample_details$s2c[, ID])
		  all_same <- function(v) {
			if(length(unique(v)) == 1) {
			  TRUE
			} else {
			  FALSE
			}
		  }
		  s2c_same <- self$sample_details$s2c %>%
			dplyr::group_by(across(all_of(ID))) %>%
		  	dplyr::summarise(across(.fns=all_same))
		  keep <- unlist(apply(s2c_same, 2, all))
		  keep[ID] <- T
		  s2c <- self$sample_details$s2c %>%
			dplyr::select(which(keep)) %>%
		  	dplyr::group_by(across(all_of(ID))) %>%
			dplyr::slice(1) %>%
			dplyr::ungroup()
		  if (!all(dim(s2c)==dim(self$sample_details$s2c))) {
			print(paste("Summing technical replicates and removing rows where",  ID, "is duplicated"))
			print(paste("Also removing columns that would be ambiguous, i.e. where summed replicates differ"))
			self$sample_details$s2c <- s2c
		  }
		}
	  	self$dge_list <- edgeR::calcNormFactors(self$dge_list) # default method is TMM
	  	# i found in at least some cases that the norm factors affect filterbyexpr
	    # seems more stringent to apply it after (less targets included)
	    # this case was one where only one sample was replicated with two libraries, perhaps related to that
	  	keep <- edgeR::filterByExpr(self$dge_list)
	   	self$dge_list <- self$dge_list[keep, ]
	   	return(invisible(self))
	},
	design = function(model) {
	  	s2c <- dplyr::filter(self$sample_details$s2c, !(sample %in% self$outliers))
	  	design <- model.matrix(model, s2c)
	  	return(design)
	},
	fit = function(model, block = NULL) {
	  	design <- self$design(model)
	  	voom_object <- limma::voom(self$dge_list, design)
	  	if (is.null(block)) {
		  	fit <- limma::lmFit(voom_object, design)
		} else {
		  	blocking_vector <- sample_details$s2c[
			  !(self$sample_details$s2c$sample %in% self$outliers),
			  colnames(self$sample_details$s2c) == block
			]
		  	cor <- self$calculate_correlation(voom_object, design, blocking_vector)
		  	fit <- limma::lmFit(
			  voom_object,
			  block = blocking_vector,
			  correlation = cor
			)
		}
		return(limma::eBayes(fit))
	},
	contrast_fit = function(fit, contrasts) {
	  	fit <- limma::contrasts.fit(fit, contrasts)
	  	return(limma::eBayes(fit))
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