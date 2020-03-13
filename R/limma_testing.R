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
	initialize = function(pathing, sample_details, reference_details) {
		self$base_path <- pathing$return_to_base()$append_to_path('quality_control')$current_path
		self$sample_details <- sample_details
		self$reference_details <- reference_details
		return(invisible(self))
	},
	tx_import = function() {
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
		colnames(txi_kallisto$counts) <- self$sample_details$s2c$sample
		return(invisible(self))
	},
	 prepare_dgelist = function(grouping_vector) {
	 	self$dge_list <- edgeR::DGEList(
		  counts = self$txi_kallisto$counts,
		  samples = sample_details$s2c,
		  group = factor(paste(
			grouping_vector,
			collapse='.')
		  )
		)
	   	keep <- edgeR::filterByExpr(self$dge_list)
	   	self$dge_list <- self$dge_list[keep, ]
	   	self$dge_list <- edgeR::calcNormFactors(self$dge_list)
	},
	fit = function(model, block = NULL) {
	  	design <- model.matrix(reformulate(model), self$dge_list$samples)
	  	voom_object <- limma::voom(self$dge_list, design)
	  	if (is.null(block)) {
		  	fit <- limma::lmfit(voom_object, design)
		} else {
		  	cor <- self$calculate_correlation(block, voom_object)
		  	fit <- limma::lmFit(voom_object, block = block, correlation = cor)
		}
		fit <- limma::eBayes(fit)
	},
	calculate_correlation =function (block, voom_object) {
		corfit <- limma::duplicateCorrelation(voom_object, design, block = block)
	  	print(paste('Estimated average correlation between block members:', corfit$cor))
	  	return(corfit$cor)
	},
	deg_results = function(fit, coeff_index, n = NULL) {
	  	# coeff_index is the indexed position of the cofficient of interest in the model,
	  	# i.e. which collumn in model.matrix
	  	if (is.null(n)) {
		  n <- length(fit$F)
		}
	  	limma::topTable(fit, coeff = coeff_index, number = n)
	},
	ndeg = function(fit) {
	  	summary(limma::decideTests(fit))
	},

))