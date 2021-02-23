# Title     : Gather sample details
# Objective : Gather sample details
# Created by: marcus
# Created on: 10/02/2020

SampleDetails <- R6::R6Class(
  classname = "SampleDetails",
  public = list(
	s2c = NULL,
	initialize = function(s2c_file_path, kallisto_dir_path, sample_quality_file_path = NULL) {
	  if (!file.exists(s2c_file_path)) {
		stop('Samples to condition matrix (s2c) not found in provided path')
	  } else if (!dir.exists(kallisto_dir_path)) {
		stop('Kallisto data directory was not found')
	  }
	  # Load sample to condition mapping
	  self$s2c <- read.table(s2c_file_path, header = TRUE, sep = ",")
	  # Add paths to kallisto ouput
      if (!('sample' %in% colnames(self$s2c))) {
        stop(paste(
          'Samples to conditions matrix must have a sample name column',
          'which should correspond to the names of the  kallisto folders'
        ))
      }
	  if ('path' %in% colnames(self$s2c)) {
		self$s2c <- dplyr::mutate(self$s2c, path = file.path(kallisto_dir_path, self$s2c[, 'path']))
	  }
	  else {
	    print('Warning, path not found in sample details, will attemp to guess it from the sample names')
        self$s2c <- dplyr::mutate(self$s2c, path = file.path(kallisto_dir_path, self$s2c[, 'sample']))
	  }
      # load sample quality details
      if (!is.null(sample_quality_file_path)) {
        if (!file.exists(sample_quality_file_path)) {
          stop('Sample quality file not found in the provided path')
        }
        sample_input <- read.table(sample_quality_file_path, sep=',', header=T)
        self$s2c$source_sample <- stringr::str_split_fixed(self$s2c$sample, '_', 2)[, 1]
        names(sample_input)[1] <- 'source_sample'
        sample_input$Date <- as.Date(sample_input$Date, format='%d.%m.%y')
        self$s2c <- dplyr::left_join(self$s2c, sample_input, by='source_sample', sort=FALSE)
        rm(sample_input)
      }
      invisible(self)
	},
    # provide factors in the order desired for default sorting
    set_factors = function(factors) {
      lapply(
        seq_along(factors),
        function(i) {
          self$s2c[, factors[i]] <- factor(
            self$s2c[, factors[i]],
            levels=unique(self$s2c[, factors[i]][order(self$s2c[, factors])])
          )
        })

      invisible(self)
    }
  )
)

ReferenceDetails <- R6::R6Class(
  classname = "ReferenceDetails",
  public = list(
    t2g = NULL,
	rRNA_genes = NULL,
  	# Transcript to gene mapping
    load_t2g = function(file_path, sep=",", header=FALSE) {
	  # Ensure input is a table with first column being transcript ID and second being gene ID
      if (!file.exists(file_path)) {
		print(file_path)
        stop('Transcripts to genes file not found in provided path')
      }
	  t2g <- read.table(file_path, sep = sep, header = FALSE)
	  if (header) {
		t2g <- t2g[-1,]
	  }
      self$t2g <- dplyr::rename(t2g, target_id = 'V1', gene_id = 'V2')
	  invisible(self)
	},
	load_rRNA_genes = function(rRNA_targets_path) {
	  rRNA_targets <- read.table(rRNA_targets_path, header=FALSE, stringsAsFactors = FALSE)
	  colnames(rRNA_targets) <- 'target_id'
	  rRNA_targets$target_id <- paste0('ref|', rRNA_targets$target_id, '|')
	  rRNA_genes <- dplyr::left_join(rRNA_targets, self$t2g, by='target_id')
	  self$rRNA_genes <- unique(rRNA_genes$gene_id)
	  invisible(self)
	}
  )
)
