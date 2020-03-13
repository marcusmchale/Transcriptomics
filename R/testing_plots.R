# Title     : Creates testing plots
# Objective : Object to gather various testing plots
# Created by: marcus
# Created on: 13/02/2020


plot_volcano_custom <- function (
	merged_results,
	file_path,
	top_n = 10,
	point_alpha = 0.2
) {
	if ('interaction_qval' %in% colnames(merged_results)) {
		res <- droplevels(dplyr::mutate(
			merged_results,
			significant = as.factor(ifelse(
				(qval.LRT <= 0.05 &	qval.WT <= 0.05),
				ifelse(
					interaction_qval > 0.05,
					'Affected (not interacting)', # Affected and no interaction with covariates
					ifelse(
						main_effect_qval <= 0.05,
						'Interacting (main effect)', # Interacting but still has main effect
						'Interacting (no main effect)' # not affected by condition
					)
				),
				'Not affected' # not affected by condition
			))
		))
		res <- res[order(
			res[, 'significant'] %in% c('Affected (not interacting)', 'Interacting (main effect)'),
			res[, 'qval.LRT'],
			decreasing=c(T,F)
		),]
		colour_map <- list(
			'Affected (not interacting)' = 'red',
			'Interacting (main effect)' = 'yellow',
			'Interacting (no main effect)' = 'blue',
			'Not affected' = 'black'
		)
		colours_vector <- as.vector(unlist(colour_map[levels(res$significant)]))

	} else {
		res <- droplevels(dplyr::mutate(
			merged_results,
			significant = as.factor(qval.LRT <= 0.05 & qval.WT <= 0.05)
		))
		res <- res[order(
			res[, 'significant'],
			res[, 'qval.LRT'],
			decreasing=c(T,F)
		),]
		colours_vector <- c('black', 'red')
	}
	if (grepl(res$target_id[1], '|')) {
		res <- dplyr::mutate(res, target_id_simple = stringr::str_extract(target_id, "(?<=\\|)(.*?)(?=\\|)"))
	} else {
		res <- dplyr::mutate(res, target_id_simple = target_id)
	}
	p <- ggplot2::ggplot(res, ggplot2::aes(b, -log10(qval.LRT), label=target_id_simple)) +
		ggplot2::geom_point(ggplot2::aes(colour = significant,), alpha = point_alpha) +
		ggrepel::geom_text_repel(
			data=res[1:top_n,],
			min.segment.length= ggplot2::unit(0,'lines')
		) +
		ggplot2::scale_colour_manual(values = colours_vector) +
		ggplot2::xlab("beta_value") +
		ggplot2::ylab("-log10(qval.LRT)") +
		ggplot2::geom_vline(xintercept = 0, colour = "black", linetype = "longdash")
	ggplot2::ggsave(file_path, plot=p)
}



plot_per_group <- function (samples, condition, results_path, grouping_property, top_n = 10, sig_value = 0.05) {
	full_results_path <- paste(batch, 'sleuth_testing', results_path, sep='/')
	out_dir <- paste(full_results_path, 'top_candidate_plots', sep='/')
	dir.create(out_dir, recursive = TRUE)
	transcript_results_files <- list.files(
		path = paste(full_results_path, 'transcript', sep='/'),
		pattern = 'results'
	)
	condition_folder <- basename(results_path)
	split_path <- function(path) {
  		if (dirname(path) %in% c("..", path)) return(basename(path))
  		return(c(basename(path), split_path(dirname(path))))
	}
	samples_folder <- split_path(dirname(results_path))[length(split_path(dirname(results_path)))]
	groups <- levels(samples[, grouping_property])
	# transcript_resilts_files and gene_results_files should be lists of same length so just indexing through them
	for (result in transcript_results_files) {
		# Get the top 10 genes and top 10 transcripts and form a consensus list of genes
		# read in the top 100 lines and only process this to get the list, should be sorted by qval.LRT already
		head_genes_results <-read.table(
			paste(full_results_path, 'gene', result, sep="/"),
			sep = '\t',
			header = TRUE,
			nrows = 100,
			stringsAsFactors = FALSE
		)
		head_transcript_results <- read.table(
			paste(full_results_path, 'transcript', result, sep="/"),
			sep = '\t',
			header = TRUE,
			nrows = 100,
			stringsAsFactors = FALSE
		)
		if ('interaction_qval' %in% colnames(head_genes_results)) {
			head_genes_results <- dplyr::filter(
				head_genes_results,
				(qval.LRT <= sig_value & qval.WT <= sig_value & interaction_qval > sig_value) |
				(qval.LRT <= sig_value & qval.WT <= sig_value & interaction_qval <= sig_value & main_effect_qval <= sig_value)
			)
			head_transcript_results <- dplyr::filter(
				head_transcript_results,
				(qval.LRT <= sig_value & qval.WT <= sig_value & interaction_qval > sig_value) |
				(qval.LRT <= sig_value & qval.WT <= sig_value & interaction_qval <= sig_value & main_effect_qval <= sig_value)
			)
		}
		head_genes <- head_genes_results$target_id[1:top_n]
		head_transcript_genes <- head_transcript_results$gene_id.LRT[1:top_n]
		genes <- unique(c(head_genes, head_transcript_genes))
		transcript_results <- dplyr::filter(
			read.table(paste(full_results_path, 'transcript', result, sep="/"), sep= '\t', header=TRUE),
			gene_id.LRT %in% genes
		)
		# aggregate transcript level data for these genes
		# from each group modelled independently
		# in testing for differences between same levels of same condition
		per_group_data <- NULL
		for (group in groups) {
			file_to_load <- paste(
				batch,
				'sleuth_testing',
				samples_folder,
				paste0('each_', grouping_property),
				group,
				condition_folder,
				'transcript',
				result,
				sep="/"
			)
			if (file.exists(file_to_load)) {
				group_temp <- dplyr::filter(
					read.table(
						file_to_load,
						header=TRUE,
						sep= '\t'
					),
					gene_id.LRT %in% genes
				)[c('target_id', 'b','se_b')]
				if (grouping_property == 'variety') {
					group_temp <- cbind(
						"group" = group,
						"pedigree" = unique(samples[samples[, grouping_property]==group,]$pedigree),
						"family" =  unique(samples[samples[, grouping_property]==group,]$family),
						group_temp
					)
				} else { # grouping property = 'family'
					group_temp <- cbind(
						"group" = group,
						"family" =  unique(samples[samples[, grouping_property]==group,]$family),
						group_temp
					)
				}
				per_group_data <- rbind(per_group_data, group_temp)
				rm(group_temp)
			}
		}
		# enforce the levels of pedigree for display purposes in facet grid (F1 in middle)
		if (grouping_property == 'variety') {
			per_group_data$pedigree <- factor(per_group_data$pedigree, levels=c("Maternal", "F1", "Paternal"))
		}
		for (gene in genes) {
			if (length(levels(per_group_data$group)) > 1) {
				transcripts <- t2g[t2g$gene_id==gene,]$target_id
				gene_group_df <- per_group_data[per_group_data$target_id %in% transcripts,]
				if (grouping_property == 'variety') {
					facets <- c('pedigree', 'family')
 				} else if (grouping_property == 'family') {
					facets <- 'family'
				}
				p <- transcripts_plot(
					data = gene_group_df,
					intercept = levels(samples[,condition])[1],
					facets = facets,
					title_text = paste(gene, 'transcripts', sep = ' ')
				)


				file_path <- paste(
					out_dir,
					sub('_results.txt', '', result),
					paste0(gene, '_per_', grouping_property , '.png', sep=""),
					sep = '/'
				)
				dir.create(dirname(file_path), recursive = TRUE)
				print(file_path)
				ggplot2::ggsave(file=file_path, plot=p)
			}
			# Plot the fold change from the all groups together model (c.f per group member facets above)
			transcripts <- t2g[t2g$gene_id==gene,]$target_id
			results_gene_df <- transcript_results[transcript_results$target_id%in%transcripts,][c('target_id','b','se_b', 'mean_obs.LRT')]
			p_all <- transcripts_plot(
				results_gene_df,
				title_text = paste(gene, 'transcripts', sep = ' ')
			)
			file_path <- paste(
				out_dir,
				sub('_results.txt', '', result),
				paste0(gene, '.png'),
				sep = '/'
			)
			dir.create(dirname(file_path), recursive = TRUE)
			ggplot2::ggsave(file=file_path, plot=p_all)
		}
	}
}


plot_per_group_new <- function (
	samples,
	condition,
	results_path,
	grouping_property,
	gene_aggregation = FALSE,
	plot_count = 10,
	sig_value = 0.05
) {
	# prepare a full path including batch and subfolder ('sleuth_testing')
	full_results_path <- paste(batch, 'sleuth_testing', results_path, sep='/')
	# get the path for the sample set assessed to retrieve per group data
	samples_folder <- split_path(dirname(results_path))[length(split_path(dirname(results_path)))]
	# the names to check for per group analysis within this sample set
	groups <- levels(samples[, grouping_property])
	# and the condition folder to find within that per group analysis
	condition_folder <- basename(results_path)
	# Now we want to select candidates to plot
	for (gt in c('gene', 'transcript')) {

		results_files <- list.files(
			path = paste(full_results_path, gt, sep='/'),
			pattern = 'results'
		)


		for (result in results_files) {
			# don't read in the whole file, just first n rows, should be sorted appropriately already
			head_results <- read.table(
				paste(full_results_path, 'transcript', result, sep="/"),
				sep = '\t',
				header = TRUE,
				nrows = 100,
				stringsAsFactors = FALSE
			)
		}
	}
	# transcript_resilts_files and gene_results_files should be lists of same length so just indexing through them
	for (result in transcript_results_files) {
		# Get the top 10 genes and top 10 transcripts and form a consensus list of genes
		# read in the top 100 lines and only process this to get the list, should be sorted by qval.LRT already
		head_genes_results <-read.table(
			paste(full_results_path, 'gene', result, sep="/"),
			sep = '\t',
			header = TRUE,
			nrows = 100,
			stringsAsFactors = FALSE
		)
		head_transcript_results <- read.table(
			paste(full_results_path, 'transcript', result, sep="/"),
			sep = '\t',
			header = TRUE,
			nrows = 100,
			stringsAsFactors = FALSE
		)
		if ('interaction_qval' %in% colnames(head_genes_results)) {
			head_genes_results <- dplyr::filter(
				head_genes_results,
				(qval.LRT <= sig_value & qval.WT <= sig_value & interaction_qval > sig_value) |
				(qval.LRT <= sig_value & qval.WT <= sig_value & interaction_qval <= sig_value & main_effect_qval <= sig_value)
			)
			head_transcript_results <- dplyr::filter(
				head_transcript_results,
				(qval.LRT <= sig_value & qval.WT <= sig_value & interaction_qval > sig_value) |
				(qval.LRT <= sig_value & qval.WT <= sig_value & interaction_qval <= sig_value & main_effect_qval <= sig_value)
			)
		}
		head_genes <- head_genes_results$target_id[1:top_n]
		head_transcript_genes <- head_transcript_results$gene_id.LRT[1:top_n]
		genes <- unique(c(head_genes, head_transcript_genes))
		transcript_results <- dplyr::filter(
			read.table(paste(full_results_path, 'transcript', result, sep="/"), sep= '\t', header=TRUE),
			gene_id.LRT %in% genes
		)
		# aggregate transcript level data for these genes
		# from each group modelled independently
		# in testing for differences between same levels of same condition
		per_group_data <- NULL
		for (group in groups) {
			file_to_load <- paste(
				batch,
				'sleuth_testing',
				samples_folder,
				paste0('each_', grouping_property),
				group,
				condition_folder,
				'transcript',
				result,
				sep="/"
			)
			if (file.exists(file_to_load)) {
				group_temp <- dplyr::filter(
					read.table(
						file_to_load,
						header=TRUE,
						sep= '\t'
					),
					gene_id.LRT %in% genes
				)[c('target_id', 'b','se_b')]
				if (grouping_property == 'variety') {
					group_temp <- cbind(
						"group" = group,
						"pedigree" = unique(samples[samples[, grouping_property]==group,]$pedigree),
						"family" =  unique(samples[samples[, grouping_property]==group,]$family),
						group_temp
					)
				} else { # grouping property = 'family'
					group_temp <- cbind(
						"group" = group,
						"family" =  unique(samples[samples[, grouping_property]==group,]$family),
						group_temp
					)
				}
				per_group_data <- rbind(per_group_data, group_temp)
				rm(group_temp)
			}
		}
		# enforce the levels of pedigree for display purposes in facet grid (F1 in middle)
		if (grouping_property == 'variety') {
			per_group_data$pedigree <- factor(per_group_data$pedigree, levels=c("Maternal", "F1", "Paternal"))
		}
		for (gene in genes) {
			if (length(levels(per_group_data$group)) > 1) {
				transcripts <- t2g[t2g$gene_id==gene,]$target_id
				gene_group_df <- per_group_data[per_group_data$target_id %in% transcripts,]
				if (grouping_property == 'variety') {
					facets <- c('pedigree', 'family')
 				} else if (grouping_property == 'family') {
					facets <- 'family'
				}
				p <- transcripts_plot(
					data = gene_group_df,
					intercept = levels(samples[,condition])[1],
					facets = facets,
					title_text = paste(gene, 'transcripts', sep = ' ')
				)


				file_path <- paste(
					out_dir,
					sub('_results.txt', '', result),
					paste0(gene, '_per_', grouping_property , '.png', sep=""),
					sep = '/'
				)
				dir.create(dirname(file_path), recursive = TRUE)
				print(file_path)
				ggplot2::ggsave(file=file_path, plot=p)
			}
			# Plot the fold change from the all groups together model (c.f per group member facets above)
			transcripts <- t2g[t2g$gene_id==gene,]$target_id
			results_gene_df <- transcript_results[transcript_results$target_id%in%transcripts,][c('target_id','b','se_b', 'mean_obs.LRT')]
			p_all <- transcripts_plot(
				results_gene_df,
				title_text = paste(gene, 'transcripts', sep = ' ')
			)
			file_path <- paste(
				out_dir,
				sub('_results.txt', '', result),
				paste0(gene, '.png'),
				sep = '/'
			)
			dir.create(dirname(file_path), recursive = TRUE)
			ggplot2::ggsave(file=file_path, plot=p_all)
		}
	}
}

transcripts_plot <- function(
	data,
	intercept,
	facets = NULL,
	transcript_axis_labels = FALSE,
	title_text = paste(datatarget_id, collapse=', '),
	x = 'target_id'
) {
	p <- ggplot2::ggplot(data, ggplot2::aes_string(x=x, y='b', fill=x)) +
		ggplot2::ggtitle(title_text) +
		ggplot2::ylab(paste("log2 fold change relative to", intercept, sep=' ')) +
		ggplot2::xlab(ggplot2::element_blank()) +
		ggplot2::theme(axis.text.x=ggplot2::element_blank()) +
		ggplot2::geom_bar(position="dodge", stat="identity") +
		ggplot2::geom_errorbar(ggplot2::aes(ymin=b-se_b, ymax=b+se_b), width=.2,position=ggplot2::position_dodge(.9)) +
		ggplot2::geom_text(ggplot2::aes(label=round(2^mean_obs.LRT)), position = ggplot2::position_stack(vjust = 0.5)) +
		ggplot2::theme(legend.title = ggplot2::element_text(paste(x)))
	if (transcript_axis_labels) {
		p <- p + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=50, hjust=1))
	} else {
		p <- p + ggplot2::theme(axis.text.x=ggplot2::element_blank())
	}
	if (length(facets) == 1) {
		p <- p + ggplot2::facet_wrap(reformulate(facets), scales = 'fixed')
	} else if (length(facets == 2)) {
		p <- p + ggplot2::facet_grid(reformulate(facets), scales = 'fixed')
	}
	return(p)
}

plot_expression_classes <- function (samples, results_path, sample_n = 10) {
	full_results_path <- paste(batch, 'sleuth_testing', results_path, sep ='/')
	condition_folder <- basename(results_path)
	if (grepl('pedigree', results_path)) {
		family_members <- c('F1', 'Maternal', 'Paternal')
		expression_classes <- expression_classes_pedigree
	} else {
		family_members <- c('0_Hybrid', '1_Inbred', '2_Wild')
		expression_classes <- expression_classes_heterozygosity
	}
	split_path <- function(path) {
  		if (dirname(path) %in% c("..", path)) return(basename(path))
  		return(c(basename(path), split_path(dirname(path))))
	}
	samples_folder <- split_path(dirname(results_path))[length(split_path(dirname(results_path)))]
	families <- levels(samples[, 'family'])
	for (gt in c('gene', 'transcript')) {
		candidate_files <- list.files(
			path = paste(full_results_path, gt, 'expression_classes', sep='/'),
			pattern = '^7_|^8_|^9_|^10_|^11_|^12_'
		)
		for (i in seq_along(candidate_files)) {
			targets <-scan(
				paste(full_results_path, gt, 'expression_classes', candidate_files[i], sep="/"),
				what=character()
			)
			targets <- sample(targets, min(sample_n, length(targets)))
			expression_class_code <- regmatches(candidate_files[i], regexpr("[^_]*", candidate_files[i]))
			expression_class_text <- expression_classes[expression_classes$code == expression_class_code,]$description
			if (length(targets) != 0) {
				per_family_data <- NULL
				if (length(families) > 1) {
					for (family in families) {
						hybrid_inbred_results <- paste(
							batch,
							'sleuth_testing',
							samples_folder,
							'each_family',
							family,
							condition_folder,
							gt,
							paste(family_members[1], 'v', family_members[2], 'results.txt', sep='_'),
							sep="/"
						)
						hybrid_wild_results <- paste(
							batch,
							'sleuth_testing',
							samples_folder,
							'each_family',
							family,
							condition_folder,
							gt,
							paste(family_members[1], 'v', family_members[3], 'results.txt', sep='_'),
							sep="/"
						)
						if (file.exists(hybrid_inbred_results)) {
							hybrid_inbred <- dplyr::filter(
								read.table(
									hybrid_inbred_results,
									header=TRUE,
									sep= '\t'
								),
								target_id %in% targets
							)[c('target_id', 'b','se_b')] %>%
								dplyr::mutate(parent = family_members[2], family = family)
						} else {
							hybrid_inbred <- NULL
						}
						if (file.exists(hybrid_wild_results)) {
							hybrid_wild <- dplyr::filter(
								read.table(
									hybrid_wild_results,
									header=TRUE,
									sep= '\t'
								),
								target_id %in% targets
							)[c('target_id', 'b','se_b')] %>%
								dplyr::mutate(parent = family_members[3], family = family)
						} else {
							hybrid_wild <- NULL
						}
						per_family_data <- rbind(per_family_data, hybrid_inbred, hybrid_wild)
						rm(hybrid_inbred, hybrid_inbred_results, hybrid_wild, hybrid_wild_results)
					}
					# parent and family as factors
					per_family_data$parent <- factor(per_family_data$parent)
					per_family_data$family <- factor(per_family_data$family)
				}
				# Also collect data for all families
				hybrid_inbred_results <- paste(
						full_results_path,
						gt,
						paste(family_members[1], 'v', family_members[2], 'results.txt', sep='_'),
						sep="/"
					)
				hybrid_wild_results <- paste(
						full_results_path,
						gt,
						paste(family_members[1], 'v', family_members[3], 'results.txt', sep='_'),
						sep="/"
					)
				if (file.exists(hybrid_inbred_results)) {
					hybrid_inbred <- dplyr::filter(read.table(
						hybrid_inbred_results,
						header= TRUE,
						sep= '\t'
					),	target_id %in% targets
					)[c('target_id', 'b','se_b')] %>%
						dplyr::mutate(parent = family_members[2])
				} else {
					hybrid_inbred <- NULL
				}
				if (file.exists(hybrid_wild_results)) {
					hybrid_wild <- dplyr::filter(read.table(
						hybrid_wild_results,
						header= TRUE,
						sep= '\t'
					),	target_id %in% targets
					)[c('target_id', 'b','se_b')] %>%
						dplyr::mutate(parent = family_members[3])
				} else {
					hybrid_wild <- NULL
				}
				all_family_data <- rbind(hybrid_inbred, hybrid_wild)
				rm(hybrid_inbred, hybrid_inbred_results, hybrid_wild, hybrid_wild_results)
				for (target in targets) {
					if (length(levels(as.factor(per_family_data$family))) > 1) {
						target_df <- per_family_data[per_family_data$target_id == target, ]

						p <- transcipts_plot(target_df, x='parent', facets = 'family', title_text = '')


						p <- ggplot2::ggplot(target_df, ggplot2::aes(x=parent, y=b, fill=parent)) +
							ggplot2::ggtitle(paste(target, '(expression class', expression_class_text, ') Hybrid Vs Parent', sep=" ")) +
							ggplot2::ylab("log2FC") +
							ggplot2::xlab(ggplot2::element_blank()) +
							ggplot2::geom_bar(position="dodge", stat="identity") +
							ggplot2::geom_errorbar(ggplot2::aes(ymin=b-se_b, ymax=b+se_b), width=.2,position=ggplot2::position_dodge(.9)) +
							ggplot2::facet_wrap(~family, scales = 'fixed')
						file_path <- paste(
							full_results_path,
							gt,
							'expression_classes',
							paste0(sub('.txt', '', candidate_files[i]), '_example_plots'),
							paste0(target, '_per_family.png', sep=""),
							sep = '/'
						)
						dir.create(dirname(file_path), recursive = TRUE)
						print(file_path)
						ggplot2::ggsave(file=file_path, plot=p)
						rm(target_df)
					}
					if (!is.null(all_family_data)) {
						all_family_data$parent <- factor(all_family_data$parent)
						# Plot the fold change from the all families together model (c.f per family facets above)
						target_df <- all_family_data[all_family_data$target_id == target,]
						p_all <-  ggplot2::ggplot(target_df, ggplot2::aes(x=parent, y=b, fill=parent)) +
							ggplot2::ggtitle(paste(target, '(expression class', expression_class_text, ') Hybrid Vs Parent', sep=" ")) +
							ggplot2::ylab("log2FC") +
							ggplot2::xlab(ggplot2::element_blank()) +
							ggplot2::geom_bar(position="dodge", stat="identity") +
							ggplot2::geom_errorbar(ggplot2::aes(ymin=b-se_b, ymax=b+se_b), width=.2,position=ggplot2::position_dodge(.9)) +
							ggplot2::theme(axis.text.x=ggplot2::element_text(angle=50, hjust=1))
						file_path <- paste(
								full_results_path,
								gt,
								'expression_classes',
								paste0(sub('.txt', '', candidate_files[i]), '_example_plots'),
								paste0(target, '.png', sep=""),
								sep = '/'
						)
						dir.create(dirname(file_path), recursive = TRUE)
						ggplot2::ggsave(file=file_path, plot=p_all)
					}
				}
			}
		}
	}
}




