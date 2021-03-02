# Title     : Quality control
# Objective : Applies quality control procedure and renders figures to examine the results
# Created by: marcus
# Created on: 10/02/2020

QC <- R6::R6Class("QC", list(
	base_path = NULL,
	sample_details = NULL,
  	reference_details = NULL,
	t2g = NULL,
	txi_kallisto = NULL,
	normalised_counts = NULL,
	outliers = NULL,
	n_common_genes = NULL,
	detection_threshold = 5,
	missing_samples_threshold = 1,
	pc_data = list(),
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
		tx2gene = self$reference_details$t2g
	  )
	  colnames(self$txi_kallisto$counts) <- self$sample_details$s2c$sample
	  # add count data to sample details table
	  self$sample_details$s2c <- dplyr::mutate(self$sample_details$s2c, est_count = colSums(self$txi_kallisto$counts))
	  write.table(
		  data.frame(self$txi_kallisto$counts) %>% tibble::rownames_to_column(var='target_id'),
		  file = file.path(pathing$base_path, 'raw_counts.tsv'),
		  row.names = FALSE,
		  quote = FALSE,
		  sep = '\t'
		)
	  return(invisible(self))
	},
	calculate_normalised_counts = function () {
	  if (is.null(self$txi_kallisto)) {
		self$tx_import()
	  }
	  s2c <- self$sample_details$s2c
	  normalised_counts <- DESeq2::DESeqDataSetFromTximport(
		  self$txi_kallisto,
		  s2c,
		  reformulate('sample')
	  ) %>%
		DESeq2::estimateSizeFactors() %>%
		BiocGenerics::counts(normalized=TRUE)
	  normalised_counts <- data.frame(normalised_counts)
	  colnames(normalised_counts) <- s2c$sample
	  self$normalised_counts <- normalised_counts
	  write.table(
		  normalised_counts %>% tibble::rownames_to_column(var='target_id'),
		  file = file.path(pathing$base_path, 'normalised_counts.tsv'),
		  row.names = FALSE,
		  quote = FALSE,
		  sep = '\t'
		)
	  return(invisible(self))
	},
	calculate_ribosomal_read_percent = function () {
	  if (is.null(self$txi_kallisto)) {
		self$sample_details$tx_import()
	  }
	  if (is.null(reference_details$rRNA_genes)) {
		stop('Import rRNA genes first (load_rRNA_genes in ReferenceDetails)')
	  }
	  rRNA_genes <- reference_details$rRNA_genes
	  # we don't use normalised counts as we are just looking at proportions anyway
	  countdata_rRNA <- self$txi_kallisto$counts[row.names(self$txi_kallisto$counts) %in% as.character(rRNA_genes),]
	  self$sample_details$s2c$rRNA_reads_percent <- 100*(colSums(countdata_rRNA)/self$sample_details$s2c$est_count)
	  return(invisible(self))
	},
	calculate_missing_genes = function () {
	  # count genes detected with more than detection threshold est_counts in
	  # as many libraries as determined by the number of libraries minus the samples_threshold (after normalisation)
	  if (is.null(self$txi_kallisto)) {
		self$sample_details$tx_import()
	  }
	  if (is.null(self$normalised_counts)) {
		self$calculate_normalised_counts()
	  }
	  s2c <- self$sample_details$s2c
	  common_genes <- self$normalised_counts[
		rowSums(self$normalised_counts >= self$detection_threshold) >= length(s2c$sample)-self$missing_samples_threshold,
	  ]
	  self$n_common_genes <- nrow(common_genes)
	  self$sample_details$s2c$missing_genes <- colSums(
		common_genes <= self$detection_threshold
	  )
	  return(invisible(self))
	},
	calculate_n_genes = function () {
	  # count genes detected with more than detection threshold est_counts
	  if (is.null(self$txi_kallisto)) {
		self$sample_details$tx_import()
	  }
	  if (is.null(self$normalised_counts)) {
		self$calculate_normalised_counts()
	  }
	  self$sample_details$s2c$n_genes <- colSums(self$normalised_counts >= self$detection_threshold)
	  return(invisible(self))
	},
	set_outliers = function(qc_parameters) {
	  calculation_functions <- list(
		'rRNA_reads_percent'=self$calculate_rRNA_reads_percent,
		'n_genes'=self$calculate_n_genes,
		'missing_genes'=self$calculate_missing_genes
	  )
	  outliers <- vector()
	  for (q in qc_parameters) {
		if (!(qc_parameters %in% colnames(self$sample_details$s2c))) {
		  calculation_functions[[q]]()
		}
		q_outliers <- self$sample_details$s2c$sample[
		  self$sample_details$s2c[[q]] %in% boxplot.stats(self$sample_details$s2c[[q]])$out
		]
		outliers <- c(
		  as.character(outliers),
		  as.character(q_outliers)
		)
	  }
	  self$outliers <- unique(outliers)
	  return(invisible(self))
	},
	# PCA analysis without current outliers, stored in list of PCA results by name
	calculate_pc = function(
	  sample_set, # a name for the sample se, used for scree plot and to access data in pc_data list
	  outliers = self$outliers,
	  n_dims = min(c(8,length(sample_details$s2c$sample)-1)),
	  scree_plot = TRUE
	) {
	  data <- DESeq2::vst(round(self$txi_kallisto$counts))
	  colnames(data) <- self$sample_details$s2c$sample
	  data <- data[, !(colnames(data) %in% outliers)]
	  s2c <- self$sample_details$s2c[!(self$sample_details$s2c$sample %in% outliers),]
	  pcDat <- prcomp(t(data), center=TRUE, scale=FALSE)
	  if (scree_plot) {
		png(paste(self$base_path, paste0('screeplot_', sample_set, '.png'), sep='/'))
		print(screeplot(pcDat))
		dev.off()
	  }
	  importance <- summary(pcDat)$importance[2,]
	  pc_out <- dplyr::mutate(as.data.frame(pcDat$x)[,1:n_dims], sample = s2c$sample) %>%
		dplyr::rename_at(
			names(importance[1:n_dims]),
			~ paste(names(importance[1:n_dims]),
			paste0(round(importance[1:n_dims]*100), '%'))
		) %>%
		tidyr::gather(key = pc, value = coordinates, 1:n_dims) %>%
		dplyr::left_join(., s2c)
	  self$pc_data[sample_set] = list(dplyr::inner_join(
		pc_out,
		pc_out,
		by=colnames(pc_out)[!colnames(pc_out) %in% c('pc','coordinates')]
	  ))
	  return(invisible(self))
	},
	get_log_counts = function (common_genes_only = FALSE) {
	  # If detection threshold is set then only genes detected with more than this frequency in
	  # as many libraries as determined by the number of libraries minus the missing samples_threshold
	  s2c <- self$sample_details$s2c
	  if (common_genes_only) {
	  	normalised_counts <- self$normalised_counts[
		  rowSums(self$normalised_counts >= self$detection_threshold) >= length(s2c$sample)-self$missing_samples_threshold,
		]
	  } else {
		normalised_counts <- self$normalised_counts
	  }
	  logcounts <- log2(normalised_counts + 1) %>%
		tibble::rownames_to_column(var='target_id') %>%
		tidyr::pivot_longer(
			cols=colnames(normalised_counts),
			names_to='sample',
			values_to='log2_est_count'
		) %>%
		dplyr::left_join(s2c)
	  return(logcounts)
	}
))


QCPlots <- R6::R6Class("QCPlots", list(
	qc = NULL,
	width = 297, # A4 long side
	height = 210, # A4 short side
	units = 'mm',
	initialize = function(qc) {
	  self$qc <- qc
	  return(invisible(self))
	},
	depth_plot = function(colour = NULL, fill = NULL, custom_colours = NULL, save_to_file = TRUE) {
	  s2c <- self$qc$sample_details$s2c
	  s2c$sample <- factor(s2c$sample, levels=s2c$sample[order(s2c[, fill], s2c[, colour])])
	  p <- ggplot2::ggplot(
		s2c,
		ggplot2::aes(x=s2c$sample)
	  ) +
		ggplot2::geom_bar(
			ggplot2::aes(
			  weight=s2c$est_count,
			  colour=s2c[, colour],
			  fill=s2c[, fill]
			)
		) +
		ggplot2::labs(colour = as.character(colour), fill = as.character(fill)) +
		ggplot2::geom_hline(yintercept=mean(s2c$est_count)) +
		ggplot2::ylab('Estimated read count') +
		ggplot2::xlab('Library') +
		ggplot2::ggtitle('Estimated read count per sequenced library') +
	  	ggplot2::theme(
			  axis.text= ggplot2::element_text(angle=45, hjust=1),
			  plot.margin = ggplot2::unit(c(20,20,20,20), 'mm')
		)
	  if (!is.null(custom_colours)){
	  		p <- p + ggplot2::scale_colour_manual(values=custom_colours)
	  }
	  if (save_to_file) {
		ggplot2::ggsave(
		  file=file.path(self$qc$base_path, 'depth.png'),
		  plot=p,
		  width=self$width,
		  height=self$height,
		  units=self$units
	  	)
		return(invisible(self))
	  } else {
	  	return(p)
	  }
	},
	rRNA_plot = function(fill = NULL, shape = NULL, custom_shapes = NULL, save_to_file = TRUE) {
	  s2c <- self$qc$sample_details$s2c
	  if (!('rRNA_reads_percent' %in% colnames(s2c))) {
		self$qc$calculate_ribosomal_read_percent()
		s2c <- self$qc$sample_details$s2c
	  }
	  s2c$sample <- factor(s2c$sample, levels=s2c$sample[order(-s2c$rRNA_reads_percent)])
	  box_stats <- boxplot.stats(s2c$rRNA_reads_percent)$stats
	  p <- ggplot2::ggplot(
			s2c,
			ggplot2::aes(
				x = s2c$sample,
				y = s2c$rRNA_reads_percent,
				fill = if(is.null(fill)) 'black' else s2c[, fill],
				shape = if(is.null(shape))'circle' else s2c[, shape]
			)
	  ) +
		ggplot2::geom_hline(yintercept=box_stats[1],linetype='dashed', color='red') +
		ggplot2::geom_hline(yintercept=box_stats[5],linetype='dashed', color='red') +
		ggplot2::geom_point(size=3) +
		#ggrepel::geom_text_repel(
		#  ggplot2::aes(
		#	label=s2c[,fill]
		#  ),
		#  min.segment.length= ggplot2::unit(0,'lines')
		#) +
		ggplot2::ggtitle('Percent of reads from rRNA genes per sequenced library') +
		ggplot2::xlab('Library') +
		ggplot2::ylab('Percentage of reads from rRNA genes') +
		ggplot2::labs(fill = as.character(fill), shape = as.character(shape), size = ggplot2::element_blank()) +
		ggplot2::guides(fill = ggplot2::guide_legend(override.aes=list(shape=22)))+
		ggplot2::theme(
		  axis.text= ggplot2::element_text(angle=45, hjust=1),
		  plot.margin = ggplot2::unit(c(20,20,20,20), 'mm')
		)
	  if (!is.null(custom_shapes)){
	  		p <- p + ggplot2::scale_shape_manual(values=custom_shapes)
	  }
	  if (save_to_file) {
		ggplot2::ggsave(
		  file=file.path(self$qc$base_path, 'rRNA.png'),
		  plot=p,
		  width=self$width,
		  height=self$height,
		  units=self$units
		)
		return(invisible(self))
	  } else {
		return(p)
	  }
	},
	missing_genes_plot = function(fill = NULL, shape = NULL, custom_shapes = NULL, save_to_file = TRUE) {
	  s2c <- self$qc$sample_details$s2c
	  if (!('missing_genes' %in% colnames(s2c))) {
		self$qc$calculate_missing_genes()
		s2c <- self$qc$sample_details$s2c
	  }
	  s2c$sample <- factor(s2c$sample, levels=s2c$sample[order(-s2c$missing_genes)])
	  box_stats <- boxplot.stats(s2c$missing_genes)$stats
	  p <- ggplot2::ggplot(
			s2c,
			ggplot2::aes(
				x=sample,
				y=missing_genes,
				fill=s2c[, fill],
				shape=s2c[, shape]
			)
	  ) +
		ggplot2::geom_hline(yintercept=box_stats[1],linetype='dashed', color='red') +
		ggplot2::geom_hline(yintercept=box_stats[5],linetype='dashed', color='red') +
		ggplot2::geom_point(size=3) +
		#ggrepel::geom_text_repel(
		#  ggplot2::aes(
		#	label=s2c[,fill]
		#  ),
		#  min.segment.length= ggplot2::unit(0,'lines')
		#) +
		ggplot2::ggtitle(paste0(
				self$qc$n_common_genes,
				' genes were detected (i.e estimated count > ',
				self$qc$detection_threshold,
				') in almost all (',
				length(s2c$sample)-self$qc$missing_samples_threshold, '/', length(s2c$sample),
				') libraries.'
			)) +
		ggplot2::xlab('Library') +
		ggplot2::ylab('Common genes not detected') +
		ggplot2::labs(fill = as.character(fill), shape = as.character(shape), size = ggplot2::element_blank()) +
		ggplot2::guides(fill = ggplot2::guide_legend(override.aes=list(shape=22)))+
		ggplot2::theme(
		  axis.text= ggplot2::element_text(angle=45, hjust=1),
		  plot.margin = ggplot2::unit(c(20,20,20,20), 'mm')
		)
	  if (!is.null(custom_shapes)){
	  		p <- p + ggplot2::scale_shape_manual(values=custom_shapes)
	  }
	  if (save_to_file) {
		ggplot2::ggsave(
		  file=file.path(qc$base_path, 'missing_genes.png'),
		  plot=p,
		  width=self$width,
		  height=self$height,
		  units=self$units
		)
		return(invisible(self))
	  } else {
		return(p)
	  }
	},
	n_genes_plot = function(fill = NULL, shape = NULL, custom_shapes = NULL, save_to_file = TRUE) {
	  s2c <- self$qc$sample_details$s2c
	  if (!('n_genes' %in% colnames(s2c))) {
		self$qc$calculate_n_genes()
		s2c <- self$qc$sample_details$s2c
	  }
	  s2c$sample <- factor(s2c$sample, levels=s2c$sample[order(s2c$n_genes)])
	  box_stats <- boxplot.stats(s2c$n_genes)$stats
	  p <- ggplot2::ggplot(
			s2c,
			ggplot2::aes(
				x=sample,
				y=n_genes,
				fill=s2c[, fill],
				shape=s2c[, shape]
			)
	  ) +
		ggplot2::geom_hline(yintercept=box_stats[1],linetype='dashed', color='red') +
		ggplot2::geom_hline(yintercept=box_stats[5],linetype='dashed', color='red') +
		ggplot2::geom_point(size=3) +
		ggplot2::ggtitle(paste0(
				'Number of genes detected (i.e estimated count > ',
				self$qc$detection_threshold,
				') in each library'
			)) +
		ggplot2::xlab('Library') +
		ggplot2::ylab('Genes detected') +
		ggplot2::labs(fill = as.character(fill), shape = as.character(shape), size = ggplot2::element_blank()) +
		ggplot2::guides(fill = ggplot2::guide_legend(override.aes=list(shape=22)))+
		ggplot2::theme(
		  axis.text= ggplot2::element_text(angle=45, hjust=1),
		  plot.margin = ggplot2::unit(c(20,20,20,20), 'mm')
		)
	  if (!is.null(custom_shapes)){
	  		p <- p + ggplot2::scale_shape_manual(values=custom_shapes)
	  }
	  if (save_to_file) {
		ggplot2::ggsave(
		  file=file.path(qc$base_path, 'n_genes.png'),
		  plot=p,
		  width=self$width,
		  height=self$height,
		  units=self$units
		)
		return(invisible(self))
	  } else {
		return(p)
	  }
	},
	counts_plot = function(
	  colour,
	  fill,
	  custom_colours = NULL,
	  common_genes_only = FALSE,
	  save_to_file = TRUE,
	  file_suffix = NULL
	) {
	  if (common_genes_only) {
		file_suffix <- paste0('_targets_below_', self$qc$detection_threshold, '_in_less_than_', self$qc$missing_samples_threshold, '_libraries')
	  }
	  logcounts <- data.frame(self$qc$get_log_counts(common_genes_only))
	  s2c <- self$qc$sample_details$s2c
	  logcounts$sample <- factor(logcounts$sample, levels=levels(s2c$sample))
	  p <- ggplot2::ggplot(
		logcounts,
		ggplot2::aes(x=sample, y=log2_est_count)
	  ) +
		ggplot2::geom_boxplot(ggplot2::aes(
		  colour = logcounts[, colour],
		  fill = logcounts[, fill]
		)) +
		ggplot2::ylab('Log2 estimated counts') +
		ggplot2::xlab('Library') +
		ggplot2::labs(colour = colour, fill = fill) +
		ggplot2::theme(axis.text= ggplot2::element_text(angle=50, hjust=1))
	  if (!is.null(custom_colours)){
	  		p <- p + ggplot2::scale_colour_manual(values=custom_colours)
	  }
	  if(save_to_file) {
	  	ggplot2::ggsave(
		  file=file.path(self$qc$base_path, paste0('counts', file_suffix, '.png')),
		  plot=p,
		  width=self$width,
		  height=self$height,
		  units=self$units
		)
		return(invisible(self))
	  } else {
		return(p)
	  }
	},
	pca_facets_plot = function(sample_sets, panel_factors, save_to_file = TRUE) {
	  for (sample_set in sample_sets) {
		pc_data <- data.frame(self$qc$pc_data[[sample_set]])
		for (panel_factor in panel_factors) {
		  p <- ggplot2::ggplot(
			  pc_data,
			  ggplot2::aes(
				  x=pc_data$coordinates.x,
				  y=pc_data$coordinates.y,
				  colour = pc_data[, panel_factor],
			  )
			) +
			ggplot2::geom_point(
				size = 1,
				alpha = 0.2,
				stroke = 1
			) +
			ggplot2::facet_grid(pc.x ~ pc.y) +
			ggplot2::labs(
				colour = panel_factor
			) +
			ggplot2::xlab(ggplot2::element_blank()) +
			ggplot2::ylab(ggplot2::element_blank())
		  if (save_to_file) {
			ggplot2::ggsave(
				file.path(self$qc$base_path, paste0('PCA_panel_', sample_set, '_', panel_factor, '.png')),
				plot=p,
				width=self$width,
				height=self$height,
				units=self$units
			)
			if (sample_set == tail(sample_sets, n=1) & panel_factor == tail(panel_factors, n=1)) {
			  return(invisible(self))
			}
		  } else {
			return(p)
		  }
		}
	  }
	},
	pca_plot = function (sample_set, pcx, pcy, panel_factors, custom_shapes = NULL, save_to_file = TRUE) {
	  pc_data <- data.frame(self$qc$pc_data[[sample_set]])
	  pc_data <-pc_data[
		grepl(paste0('^', pcx, ' '), pc_data$pc.x) & grepl(paste0('^', pcy, ' '), pc_data$pc.y),
	  ]
	  x_label <- unique(pc_data$pc.x)
	  y_label <- unique(pc_data$pc.y)
	  if (length(x_label) > 1 | length(y_label) > 1) {
		stop('Must provide a dataframe with only one level for pc.x and one level for pc.y')
	  }
	  p <- ggplot2::ggplot(
		pc_data,
		ggplot2::aes(
		  x=pc_data$coordinates.x,
		  y=pc_data$coordinates.y,
		  colour = if(panel_factors[1] %in% colnames(pc_data)) pc_data[, panel_factors[1]] else FALSE,
		  fill = if(panel_factors[2] %in% colnames(pc_data)) pc_data[, panel_factors[2]] else FALSE,
		  shape = if(panel_factors[3] %in% colnames(pc_data)) pc_data[, panel_factors[3]] else FALSE
		)
	  ) +
		  ggplot2::geom_point(
			  size = 3,
			  alpha = 0.5,
			  stroke = 1
		  ) +
		  ggrepel::geom_text_repel(
			ggplot2::aes(
				label=paste0(pc_data$sample)
			),
			min.segment.length= ggplot2::unit(0,'lines')
		  ) +
		  ggplot2::labs(
			  colour = panel_factors[1],
			  fill = panel_factors[2],
			  shape = panel_factors[3]
		  ) +
		  ggplot2::guides(
			  fill = if(is.na(panel_factors[1])) FALSE else ggplot2::guide_legend(override.aes=list(shape=22)),
			  colour = if(is.na(panel_factors[2])) FALSE else ggplot2::guide_legend(),
			  shape = if(is.na(panel_factors[3])) FALSE else ggplot2::guide_legend()
		  ) +
		  ggplot2::xlab(x_label) +
		  ggplot2::ylab(y_label) +
		  ggplot2::scale_fill_manual(values=c('blue','yellow','red'))
	  if (!is.null(custom_shapes)){
	  		p <- p + ggplot2::scale_shape_manual(values=custom_shapes)
	  }
	  if (save_to_file) {
		ggplot2::ggsave(
			file.path(self$qc$base_path, paste0('PCA_', sample_set, '_', pcx, '_v_', pcy, '.png')),
			plot=p,
			width=self$width,
			height=self$height,
			units=self$units
		)
		return(invisible(self))
	  } else {
		return(p)
	  }
	},
	sample_plot =function(x,  xlab, y, ylab, file_suffix = NULL, save_to_file = TRUE) {
	  s2c <- self$qc$sample_details$s2c
	  s2c <- s2c[
		!(is.na(s2c[, x])) & !(is.na(s2c[, y])),
	  ]
	  p <- function() {
		plot(s2c[,x], s2c[,y], xlab = xlab, ylab= ylab, col=c(1,1,1,1))
		text(s2c[,x],s2c[,y],labels=s2c$source_sample, col=ifelse(s2c$sample %in% self$qc$outliers, rgb(1,0,0,1), rgb(0,0,1,1)))
	  }
	  if (save_to_file) {
	  	png(file.path(self$qc$base_path, paste0(xlab, '_v_', ylab, file_suffix, '.png') ))
		p()
	  	dev.off()
		return(invisible(self))
	  } else {
		return(p())
	  }
	},
	table_plot = function(x,  xlab, y, ylab, file_suffix = NULL, save_to_file = TRUE, save_table_text = TRUE) {
	  s2c <- self$qc$sample_details$s2c
	  s2c <- s2c[
		!(is.na(s2c[, x])) & !(is.na(s2c[, y])),
	  ]
	  if (save_table_text) {
		capture.output(table(s2c[,x], s2c[,y]), file=file.path(self$qc$base_path,  paste0('table_', xlab, '_v_', ylab, file_suffix, '.txt')))
	  }
	  p <- function() {
		plot(table(s2c[,x], s2c[,y]))
	  }
	  if (save_to_file) {
	  	png(file.path(self$qc$base_path, paste0('table_plot_', xlab, '_v_', ylab, file_suffix, '.png') ))
		p()
	  	dev.off()
		return(invisible(self))
	  } else {
		return(p())
	  }
	},
	capture_anova = function(dependent_var, independent_vars_string, save_to_file = TRUE) {
	  output <- anova(lm(reformulate(independent_vars_string, dependent_var), self$qc$sample_details$s2c))
	  if(save_to_file) {
		capture.output(
		  output,
		  file = file.path(self$qc$base_path, paste0('anova_', dependent_var, '.txt'))
		)
		return(invisible(self))
	  } else {
		return(output)
	  }
	}
))