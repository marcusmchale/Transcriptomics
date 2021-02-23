# Title     : OO approach to Sleuth analysis
# Objective : Rewriting the sleuth analysis scripts into an R6 framework
# Created by: marcus
# Created on: 07/02/2020


SleuthModels <- R6::R6Class("SleuthModels", list(
	s2c = NULL,
	condition = NULL,
	covariates = NULL,
	interacting_covariates = NULL,
	reduced = NULL,
	additive = NULL,
	interaction = NULL,
	interaction_less_main = NULL,
	initialize = function(s2c, condition, covariates = NULL, interacting_covariates = NULL) {
	  self$s2c <- droplevels(s2c) # make a local copy to which we add dummy coding for contr.sum with and without main effect
	  check_factors <- function (labels) {
		for (label in labels) {
		  if (!(label %in% colnames(s2c))) {
			stop(paste(label, 'not in sample details'))
		  }
		  if (!(is.factor(s2c[, label]) | is.numeric(s2c[, label]))) {
			stop(paste(label, 'is neither a factor nor numeric in sample details'))
		  }
		  print(levels(s2c[, label]))
		}
	  }
	  check_factors(c(condition, covariates, interacting_covariates))
	  # Check if proposed interaction can actually be examined
	  # (i.e. are samples present at all levels of the interaction)
	  if (!(is.null(interacting_covariates))) {
			if (sum(table(self$s2c[, c(condition, interacting_covariates)]) == 0) != 0) {
			  table(self$s2c[, c(condition, interacting_covariates)])
			  stop(paste(
				  'Cannot handle the proposed interacting covariates for this sample set',
				  'as samples are not found at all levels of the interaction'
			  ))
		}
	  }
	  if (length(condition) > 1) {
		stop('This function only develops models to test for a single condition')
	  }
	  self$condition <- condition
	  self$covariates <- unique(c(covariates, interacting_covariates))
	  self$interacting_covariates <- interacting_covariates
	  return(invisible(self))
	},
	create_reduced = function() {
	  if (is.null(self$covariates)) {
		self$reduced <- reformulate('1')
	  } else {
		self$reduced <- reformulate(paste(self$covariates, collapse = (' + ')))
	  }
	  print("Covariates (or null) additive model (reduced): ")
	  print(model.matrix(self$reduced, self$s2c))
	  return(invisible(self))
	},
	create_additive = function () {
	  self$additive <- reformulate(paste(c(self$covariates, self$condition), collapse = (' + ')))
	  print("Covariates additive condition additive model (additive) : ")
	  print(self$additive)
	  print(model.matrix(self$additive, self$s2c))
	  return(invisible(self))
	},
	create_interaction = function () {
	  non_interacting_covariates <- self$covariates[!(self$covariates %in% self$interacting_covariates)]
	  if (length(non_interacting_covariates) > 0) {
		  interaction_string <- paste(paste(non_interacting_covariates, collapse = '+' ), ' + ')
	  } else {
		  interaction_string <- NULL
	  }
	  for (cov in self$interacting_covariates) {
		cov_x_cond <- paste0('(', cov, ' * ', self$condition, ')')
		interaction_string <- paste(c(interaction_string, cov_x_cond), collapse = ' + ')
	  }
	  self$interaction <- reformulate(interaction_string)
	  print("Covariates additive condition interaction model (interaction)")
	  print(self$interaction)
	  print(model.matrix(self$interaction, self$s2c))
	  return(invisible(self))
	},
	create_interaction_less_main = function () {
	  interaction_less_main <- NULL
	  for (covariate in self$covariates) {
		  covariate.numeric <- sapply(
			self$s2c[, covariate],
			function(i) contr.sum(length(levels(self$s2c[, covariate])))[i,]
		  )
		  if (length(levels(self$s2c[, covariate])) > 2) {
			  covariate.numeric <- t(covariate.numeric)
			  colnames(covariate.numeric) <- paste0(covariate, 1:dim(covariate.numeric)[2])
			  self$s2c <- data.frame(self$s2c, covariate.numeric)
			  if (covariate %in% self$interacting_covariates) {
				  formula_string <- paste0(
					  paste(colnames(covariate.numeric), collapse = ' + '),
					  ' + ', self$condition, ':',
					  paste(colnames(covariate.numeric), collapse = paste0(' + ', self$condition, ':'))
				  )
			  }else {
				  formula_string <- paste0(
					  paste(colnames(covariate.numeric), collapse = ' + ')
				  )
			  }
		  } else {
			  self$s2c[, paste0(covariate, '1')] <- covariate.numeric
			  if (covariate %in% self$interacting_covariates) {
				  formula_string <- paste0(
					  paste0(covariate, '1'), ' + ', paste(self$condition, paste0(covariate, '1'), sep=':')
				  )
			  } else {
				  formula_string <- paste0(covariate, '1')
			  }
		  }
		  interaction_less_main <- paste(c(interaction_less_main, formula_string), collapse = ' + ')
	  }
	  self$interaction_less_main <- reformulate(interaction_less_main)
	  print("Covariates additive condition interaction model without main effect of condition (interaction_less_main) ")
	  print(self$interaction_less_main)
	  print(model.matrix(self$interaction_less_main, self$s2c))
	  return(invisible(self))
	},
	prepare_models = function() {
	  self$create_reduced()$
	  	create_additive()
	  if (!is.null(self$interacting_covariates)) {
	  	self$create_interaction()$
	  	create_interaction_less_main()
	  }
	  return(invisible(self))
	}
))

SleuthGeneComp <- R6::R6Class("SleuthGeneComp", list(
	s2c = NULL,
	t2g = NULL,
	additive = NULL,
	interaction = NULL,
	sleuth_object = NULL,
	num_cores = NULL,
	initialize = function (s2c, t2g, additive, interaction, num_cores=NULL) {
	  if (!exists('num_cores')) {
	  	self$num_cores <- parallel::detectCores()
	  } else {
		self$num_cores <- num_cores
	  }
	  print(paste(self$num_cores, "cores"))
	  self$s2c <- s2c
	  self$t2g <- t2g
	  self$additive <- additive
	  self$interaction <- interaction
	  return(invisible(self))
	},
	fit_models = function(min_prop = 0.49) {
	  self$sleuth_object <- sleuth::sleuth_prep(
			self$s2c,
			target_mapping = self$t2g,
			gene_mode = TRUE,
			aggregation_column = 'gene_id',
			transformation_function = function(x) log2(x + 0.5),
			filter_fun = function(x) sleuth::basic_filter(row=x, min_prop = min_prop),
			num_cores = self$num_cores
	  ) %>% sleuth::sleuth_fit(self$additive, 'additive')
	  if (!(is.null(self$interaction))) {
		self$sleuth_object <- sleuth::sleuth_fit(
		  self$sleuth_object,
		  self$interaction,
		  'interaction'
		)
	  }
	  return(invisible(self))
	}
))

SleuthTranscriptComp <- R6::R6Class("SleuthTranscriptComp", list(
	s2c = NULL,
	t2g = NULL,
	reduced = NULL,
	additive = NULL,
	interaction = NULL,
	interaction_less_main = NULL,
	sleuth_object = NULL,
	num_cores = NULL,
	initialize = function (s2c, t2g, reduced, additive, interaction, interaction_less_main, num_cores=NULL) {
	  if (!exists('num_cores')) {
	  	self$num_cores <- parallel::detectCores()
	  } else {
		self$num_cores <- num_cores
	  }
	  print(paste(self$num_cores, "cores"))
	  self$s2c <- s2c
	  self$t2g <- t2g
	  self$reduced <- reduced
	  self$additive <- additive
	  self$interaction <- interaction
	  self$interaction_less_main <- interaction_less_main
	  return(invisible(self))
	},
	fit_models = function(min_prop = 0.49) {
	  self$sleuth_object <- sleuth::sleuth_prep(
			self$s2c,
			target_mapping = self$t2g,
			aggregation_column = if (is.null(self$t2g)) Null else 'gene_id',
			transformation_function = function(x) log2(x + 0.5),
			filter_fun = function(x) sleuth::basic_filter(row=x, min_prop = min_prop),
			num_cores = self$num_cores
	  ) %>%
		sleuth::sleuth_fit(self$reduced, 'reduced') %>%
		sleuth::sleuth_fit(self$additive, 'additive')
	  if (!(is.null(self$interaction))) {
		self$sleuth_object <- sleuth::sleuth_fit(self$sleuth_object, self$interaction, 'interaction') %>%
		  sleuth::sleuth_fit(self$interaction_less_main, 'interaction_less_main')
	  }
	  return(invisible(self))
	},
	perform_lrt_tests = function() {
	  self$sleuth_object <-  sleuth::sleuth_lrt(self$sleuth_object, 'reduced', 'additive')
	  if (!(is.null(self$interaction))) {
		self$sleuth_object <- sleuth::sleuth_lrt(self$sleuth_object, 'reduced', 'interaction') %>%
		  sleuth::sleuth_lrt('additive', 'interaction') %>%
		  sleuth::sleuth_lrt('interaction_less_main', 'interaction')
	  }
	  return(invisible(self))
	},
	get_lrt_tables = function() {
	  	# Likelihood ratio tests
		transcript_additive_lrt_table <- sleuth::sleuth_results(
			self$sleuth_object,
			test='reduced:additive',
			test_type='lrt',
			pval_aggregate = FALSE
		)
		if (!(is.null(self$interaction))) {
			transcript_additive_v_interaction_lrt_qval <- sleuth::sleuth_results(
				self$sleuth_object,
				test='additive:interaction',
				test_type='lrt',
				pval_aggregate = FALSE
			)[,c('target_id', 'pval', 'qval')]
			colnames(transcript_additive_v_interaction_lrt_qval) <- c('target_id', 'interaction_pval', 'interaction_qval')
			transcript_interaction_less_main_v_interaction_lrt_qval <- sleuth::sleuth_results(
				self$sleuth_object,
				test='interaction_less_main:interaction',
				test_type = 'lrt',
				pval_aggregate = FALSE
			)[,c('target_id', 'pval', 'qval')]
			colnames(transcript_interaction_less_main_v_interaction_lrt_qval) <- c('target_id', 'main_effect_pval', 'main_effect_qval')
			transcript_additive_lrt_table <- merge(
				transcript_additive_lrt_table,
				transcript_additive_v_interaction_lrt_qval,
				by='target_id',
				sort = FALSE
			)
			transcript_additive_lrt_table <- merge(
				transcript_additive_lrt_table,
				transcript_interaction_less_main_v_interaction_lrt_qval,
				by='target_id',
				sort = FALSE
			)
		}
		# We use pval aggregation from same test for the gene list
		gene_additive_lrt_table <- sleuth::sleuth_results(
			self$sleuth_object,
			test='reduced:additive',
			test_type='lrt',
			pval_aggregate = TRUE
		)
		if (!(is.null(self$interaction))) {
			gene_additive_v_interaction_lrt_qval <- sleuth::sleuth_results(
				self$sleuth_object,
				test='additive:interaction',
				test_type='lrt',
				pval_aggregate = TRUE
			)[,c('target_id', 'pval', 'qval')]
			colnames(gene_additive_v_interaction_lrt_qval) <- c('target_id', 'interaction_pval', 'interaction_qval')
		  	gene_interaction_less_main_v_interaction_lrt_qval <- sleuth::sleuth_results(
				self$sleuth_object,
				test='interaction_less_main:interaction',
				test_type = 'lrt',
				pval_aggregate = TRUE
			)[,c('target_id', 'pval', 'qval')]
			colnames(gene_interaction_less_main_v_interaction_lrt_qval) <- c('target_id', 'main_effect_pval', 'main_effect_qval')
			gene_additive_lrt_table <- merge(
				gene_additive_lrt_table,
				gene_additive_v_interaction_lrt_qval,
				by='target_id',
				sort = FALSE
			)
			gene_additive_lrt_table <- merge(
				gene_additive_lrt_table,
				gene_interaction_less_main_v_interaction_lrt_qval,
				by='target_id',
				sort = FALSE
			)
		}
	  return(list(
	  	transcript = transcript_additive_lrt_table,
		gene = gene_additive_lrt_table
	  ))
	}
))

SleuthAssessCovariates <- R6::R6Class('SleuthAssessCovariates', list(
	so = NULL,
	s2c = NULL,
	t2g = NULL,
	num_cores = NULL,
	initialize = function(s2c, t2g, outliers = NULL, num_cores=NULL) {
	  if (!exists('num_cores')) {
	  	self$num_cores <- parallel::detectCores()
	  } else {
		self$num_cores <- num_cores
	  }
	  print(paste(self$num_cores, "cores"))
	  self$s2c <- s2c
	  if (!is.null(outliers)) {
		self$s2c <- self$s2c[!(self$s2c$sample %in% outliers),]
	  }
	  self$t2g <- t2g
	  self$so <- sleuth::sleuth_prep(
		s2c,
		target_mapping = t2g,
		aggregation_column = if (is.null(t2g)) Null else 'gene_id',
		transformation_function = function(x) log2(x + 0.5),
		num_cores = self$num_cores
	  )
	  return(invisible(self))
	},
	targets_affected = function(reduced, full) {
	  	so <- self$so
		so <- sleuth::sleuth_fit(
		  	so,
		  	reformulate(reduced),
		  	'reduced'
		)
	  	so <- sleuth::sleuth_fit(
		  	so,
		  	reformulate(full),
		  	'full'
		)
	  	so <- sleuth::sleuth_lrt(
			so,
			'reduced',
			'full'
		)
	  	significant_targets <- sleuth::sleuth_results(so, test='reduced:full', test_type = 'lrt', pval_aggregate = FALSE)
		dplyr::filter(significant_targets, qval <= 0.05)$target_id
	}
))


SleuthMultComp <- R6::R6Class("SleuthMultComp", list(
	sig_value = 0.05,
	base_path = NULL,
	s2c = NULL,
	t2g = NULL,
	condition = NULL,
	initialize  = function(
	  pathing,
	  sample_details,
	  reference_details,
	  outliers = NULL,
	  extra_path = NULL
	) {
	  self$base_path <- pathing$return_to_base()$append_to_path('sleuth_testing')$current_path
	  if (!is.null(extra_path)) {
		self$base_path <- file.path(self$base_path, extra_path)
	  }
	  dir.create(self$base_path, recursive = TRUE, showWarnings = FALSE)
	  self$s2c <- sample_details$s2c[, colSums(is.na(sample_details$s2c)) == 0] # get rid of any columns with missing values
	  if (!is.null(outliers)) {
		self$s2c <- self$s2c[!(self$s2c$sample %in% outliers), ]
	  }
	  self$t2g <- reference_details$t2g
	  return(invisible(self))
	},
	calculate_min_prop = function() {
	  # Prepare a condition only design matrix to calculate the min_prop to use in filtering
	  min_s <- min(colSums(model.matrix(reformulate(self$condition), self$s2c)))
	  min_s <- ifelse(min_s == 1, 2, min_s) # ensure minimum of 2
	  min_prop <- min_s / nrow(self$s2c) # calculate as proportion of all
	  min_prop <- ifelse(min_prop >=0.3, min_prop, 0.3) # hard limit of >=30% of libraries
	  return(min_prop)
	},
	get_intercepts = function() {
	  # We need to create a set of models for each intercept for WT to give fold-change between cases
	  # https://github.com/pachterlab/sleuth/issues/156
	  cond_levels <- levels(self$s2c[ ,self$condition])
	  cond_pairs <- combn(cond_levels, 2, simplify=FALSE)
	  if (length(cond_levels) < 2) {
		print(paste('Unique values for condition:', length(cond_levels)))
		print(cond_pairs)
		stop('Need at least two levels of the condition in the input data')
	  }
	  intercepts <- cond_levels[seq_len(length(cond_levels)-1)]
	  return(intercepts)
	},
	get_wt_tables = function (sleuth_object, condition_pairs, gene = TRUE) {
	  # Prepare a list to store Wald test results from each condition pair
	  wt_results <- vector(mode="list", length=length(condition_pairs))
	  names(wt_results) <- lapply(condition_pairs, function(x) paste(x, collapse=' v '))
	  for (j in seq_along(condition_pairs)) {
		  pair <- condition_pairs[[j]]
		  pair_name <- paste(pair, collapse = ' v ')
		  print(pair_name)
		  contrast <- paste0(self$condition, pair[[2]])
		  # Gene level Wald test for condition response to extract fold change from aggregated counts
		  sleuth_object <- sleuth::sleuth_wt(sleuth_object, contrast, 'additive')
		  if (gene) {
			wt_results[[pair_name]] <- sleuth::sleuth_results(  # gene_additive_wt_table
				sleuth_object,
				test=contrast,
				test_type = 'wt',
				which_model = 'additive',
				show_all= FALSE
			)
		  } else (
			wt_results[[pair_name]] <- sleuth::sleuth_results(  # gene_additive_wt_table
				sleuth_object,
				test=contrast,
				test_type = 'wt',
				which_model = 'additive',
				show_all= FALSE,
				pval_aggregate = FALSE
			)
		  )
	  }
	  return(wt_results)
	},
	min_fc_result =  function(result, min_fc=NULL) {
	  if (is.null(min_fc)) {
		return(result)
	  }
	  print('Adjusting the qval.LRT and qval.WT to reflect the arbitrary fold-change threshold')
	  print('This is done because the filtering means less comparisons are being made')
	  min_b <- log2(min_fc)
	  # replace the qvals with BH corrected qvals having filtered the targets
	  pval.LRT <-  ifelse(abs(result$b) > min_b, result$pval.LRT, 1)
	  pval.WT <- ifelse(abs(result$b) > min_b, result$pval.WT, 1)
	  result$qval.LRT <- p.adjust(pval.LRT)
	  result$qval.WT <- p.adjust(pval.WT)
	  if ('interaction_qval' %in% colnames(result)) {
		print('Also adjusting the interaction_qval and main_effect_qval for the same reason')
		interaction_pval <- ifelse(abs(result$b) > min_b, result$interaction_pval, 1)
		main_effect_pval <- ifelse(abs(result$b) > min_b, result$main_effect_pval, 1)
		result$interaction_qval <- p.adjust(interaction_pval)
		result$interaction_qval <- p.adjust(main_effect_pval)
	  }
	  return(result)
	},
	analyse_condition = function(
	  condition,
	  covariates = NULL,
	  interacting_covariates = NULL,
	  intercepts = NULL,
	  min_fc = NULL,
	  top_n=10,
	  num_cores = NULL
   	) {
	  if (is.null(num_cores)) {
	  	num_cores <- parallel::detectCores()
	  }
	  outpath <- file.path(self$base_path, condition)
	  dir.create(outpath, showWarnings = FALSE)
	  sig_value <- self$sig_value
	  self$condition <- condition
	  sleuth_models <- SleuthModels$new(self$s2c, condition, covariates, interacting_covariates)$prepare_models()
	  s2c <- sleuth_models$s2c
	  min_prop <- self$calculate_min_prop()
	  if (is.null(intercepts)) {
	 	intercepts <- self$get_intercepts()
	  }
	  for (intercept in intercepts) {
		print(paste0('Intercept is: ', intercept))
		s2c[,condition] <- relevel(s2c[,condition], ref=intercept)
		cond_levels <- levels(s2c[, condition])
	  	condition_pairs <- combn(cond_levels, 2, simplify=FALSE)
	  	condition_pairs <- Filter(function(x) x[[1]] == intercept, condition_pairs)
		wt_tables <- vector(mode='list', length = 2)
		names(wt_tables) <- c('gene', 'transcript')
	  	gene_comp <- SleuthGeneComp$new(
				s2c,
				self$t2g,
				sleuth_models$additive,
				sleuth_models$interaction,
				num_cores = num_cores
			  )$
				fit_models(min_prop)
		wt_tables[['gene']] <- self$get_wt_tables(
			  gene_comp$sleuth_object,
			  condition_pairs
		  )
		rm(gene_comp)
		transcript_comp <- SleuthTranscriptComp$new(
			s2c,
			self$t2g,
			sleuth_models$reduced,
			sleuth_models$additive,
			sleuth_models$interaction,
			sleuth_models$interaction_less_main,
			num_cores = num_cores
		)$
		  fit_models(min_prop)$
		  perform_lrt_tests()
		lrt_tables <- transcript_comp$get_lrt_tables()
		wt_tables[['transcript']] <- self$get_wt_tables(
			transcript_comp$sleuth_object,
			condition_pairs,
			gene = FALSE
		)
		rm(transcript_comp)
		for (gt in c('gene', 'transcript')) {
		  for (pair in condition_pairs) {
		  	pair_name <- paste(pair, collapse = ' v ')
			result <- merge(
			  lrt_tables[[gt]],
			  wt_tables[[gt]][[pair_name]],
			  by = 'target_id',
			  sort = FALSE,
			  suffixes = c(".LRT", ".WT")
			)
			result <- self$min_fc_result(result, min_fc)
			if (is.null(self$interaction)) {
			  	result %>% dplyr::arrange(
					qval.LRT > sig_value & qval.WT > sig_value
			  	)
			} else {
				result %>% dplyr::arrange(
					!(interaction_qval > sig_value | (interaction_qval <= sig_value & main_effect_qval <= sig_value)),
					qval.LRT
			  	)
			}
			dir.create(file.path(outpath, gt),  showWarnings = FALSE)
		  	write.table(
			  result,
			  file = file.path(outpath, gt, paste0(paste(pair, collapse='_v_'), '_result.txt')),
			  row.names = FALSE,
			  quote = FALSE,
			  sep = '\t'
			)
			plot_volcano_custom(
				result,
				file.path(outpath, gt, paste0(paste(pair, collapse='_v_'), '_volcano.png')),
				min_fc = min_fc,
				top_n = top_n
			)
			if (is.null(self$interaction)) {
			  significant_targets <- dplyr::filter(
				result,
				qval.LRT <= 0.05 & qval.WT <=0.05
			  )$target_id
			  print(
				paste0(
				  length(significant_targets),
				  ' ',
				  gt,
				  's were affected by the condition (',
				  condition,
				  ') from level ',
				  pair[[1]], ' to ', pair[[2]], '..'
				)
			  )
			  write.table(
				  data.frame(significant_targets),
				  file = file.path(outpath, gt, paste0(paste(pair, collapse='_v_'), '_diff.txt')),
				  row.names=FALSE,
				  col.names=FALSE,
				  quote=FALSE,
				  sep="\t"
			  )
			} else {
			  # Now select significant targets
			  # These should be first selected from the additive model then tested for interactions and a main effect
			  # then included only if no interaction or interaction + main effect
			  additive_targets <- dplyr::filter(result, qval.LRT <= 0.05 & qval.WT <=0.05)$target_id
			  interaction_non_sig_targets <- dplyr::filter(result, interaction_qval > 0.05)$target_id
			  interaction_sig_targets_main_effect <- dplyr::filter(result, interaction_qval <= 0.05, main_effect_qval <= 0.05)$target_id
			  additive_targets_non_interacting <- intersect(additive_targets, interaction_non_sig_targets)
			  additive_targets_interacting_main_effect <- intersect(additive_targets, interaction_sig_targets_main_effect)
			  significant_targets <- c(additive_targets_non_interacting, additive_targets_interacting_main_effect)
			  print(
				  paste0(
					  length(additive_targets_non_interacting),
					  ' ',
					  gt,
					  's were affected by the condition (',
					  condition,
					  ') from level ',
					  paste(pair, collapse = ' to '),
					  ' and this effect did not appear to interact with covariates (',
					  paste(self$covariates, collapse = ' and '),
					  '). Another ',
					  length(additive_targets_interacting_main_effect),
					  ' ',
					  gt,
					  's responded differently at different levels of some covariates (',
					  paste(self$interacting_covariates, collapse = 'and '),
					  ') but still had a detectable main effect of the condition
					  so were included in the full "diff" set (',
					  length(significant_targets),
					  ').'
				  )
			  )
			  write.table(
				  data.frame(significant_targets),
				  file = file.path(outpath, gt, paste0(paste(pair, collapse='_v_'), '_diff.txt')),
				  row.names=FALSE,
				  col.names=FALSE,
				  quote=FALSE,
				  sep="\t"
			  )
			  write.table(
				  data.frame(additive_targets_non_interacting),
				  file = file.path(outpath, gt, paste0(paste(pair, collapse='_v_'), '_diff_non_interacting.txt')),
				  row.names=FALSE,
				  col.names=FALSE,
				  quote=FALSE,
				  sep="\t"
			  )
			  write.table(
				  data.frame(additive_targets_interacting_main_effect),
				  file = file.path(outpath, gt, paste0(paste(pair, collapse='_v_'), '_diff_interacting_but_main_effect.txt')),
				  row.names=FALSE,
				  col.names=FALSE,
				  quote=FALSE,
				  sep="\t"
			  )
			}
		  }
	  	}
	  }
	  return(invisible(self))
	}

))

