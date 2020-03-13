# Title     : TODO
# Objective : TODO
# Created by: marcus
# Created on: 13/02/2020

expression_classes_pedigree <- data.frame(
	"code"=1:12,
	"description"=c(
		'M>F>P',
		'P>F>M',
		'P>M=F',
		'F=M>P',
		'M>F=P',
		'F=P>M',
		'M>P>F',
		'P>M>F',
		'P=M>F',
		'F>M>P',
		'F>P>M',
		'F>M=P'
	),
	"file_annotation"=c(
		'MFP',
		'PFM',
		'PM_F',
		'F_MP',
		'MF_P',
		'F_PM',
		'MPF',
		'PMF',
		'P_MF',
		'FMP',
		'FPM',
		'FM_P'
	)
)

expression_classes_heterozygosity <- data.frame(
	"code"=1:12,
	"description"=c(
		'I>H>W',
		'W>H>I',
		'W>I=H',
		'H=I>W',
		'I>H=W',
		'H=W>I',
		'I>W>H',
		'W>I>H',
		'W=I>H',
		'H>I>W',
		'H>W>I',
		'H>I=W'
	),
	"file_annotation"=c(
		'IHW',
		'WHI',
		'WI_H',
		'H_IW',
		'IH_W',
		'H_WI',
		'IWH',
		'WIH',
		'W_IH',
		'HIW',
		'HWI',
		'HI_W'
	)
)

classify_expression <- function(
	results_path,
	parent_labels=NULL,
	include_main_effects = FALSE
) {
	results_path <- paste(batch, 'sleuth_testing', results_path,sep='/')
	print(results_path)
	if (is.null(parent_labels)) {
		if (grepl('pedigree', results_path)) {
			parent_labels <- c('F1', 'Maternal', 'Paternal')
			expression_classes <- expression_classes_pedigree
			condition <- 'pedigree'
		} else if (grepl('heterozygosity', results_path)) {
			parent_labels <- c('0_Hybrid', '1_Inbred', '2_Wild')
			expression_classes <- expression_classes_heterozygosity
			condition <- 'heterozygosity'
		} else {
			print (
				'Parent labels not identified from input folder and need to be set (list with 3 items, in alphabetical order)'
			)
			return()
		}
	}
	print(parent_labels)
	for (g in c('transcript', 'gene')) {
		FvM_path <- paste(
			results_path,
			g,
			paste0(parent_labels[[1]], '_v_', parent_labels[[2]], '_results.txt'),
			sep='/'
		)
		FvP_path <- paste(
			results_path,
			g,
			paste0(parent_labels[[1]], '_v_', parent_labels[[3]], '_results.txt'),
			sep='/'
		)
		MvP_path <- paste(
			results_path,
			g,
			paste0(parent_labels[[2]], '_v_', parent_labels[[3]], '_results.txt'),
			sep='/'
		)
		if (!(file.exists(FvM_path)) | !(file.exists(FvP_path)) | !(file.exists(MvP_path))) {
			print(paste(
				'Cannot classify expression into heterotic classes as some levels of',
				condition,
				'are not found in the sample set.'
			))
			return('Classification failed')
		}
		FvM <- dplyr::filter(
			read.table(
				FvM_path, header = TRUE
			), qval.LRT <= 0.05
		)
		FvP <-	dplyr::filter(
			read.table(
				FvP_path, header = TRUE
			), qval.LRT <= 0.05
		)
		MvP <- dplyr::filter(
			read.table(
				MvP_path, header = TRUE
			), qval.LRT <= 0.05
		)
		if (!('interaction_qval' %in% colnames(FvM))) {
			print(paste(
					'Caution: No interaction effect between families and',
					condition,
					'could be evaluated for this sample set.',
					'As such many of the expression class candidates are likely to be influenced by variation',
					'that is specific to a single family'
			))
		} else if(include_main_effects) {
			FvM <- dplyr::filter(FvM, interaction_qval > 0.05 | (interaction_qval <=0.05 & main_effect_qval <= 0.05))
			FvP <- dplyr::filter(FvP, interaction_qval > 0.05 | (interaction_qval <=0.05 & main_effect_qval <= 0.05))
			MvP <- dplyr::filter(MvP, interaction_qval > 0.05 | (interaction_qval <=0.05 & main_effect_qval <= 0.05))
		} else {
			FvM <- dplyr::filter(FvM, interaction_qval > 0.05)
			FvP <- dplyr::filter(FvP, interaction_qval > 0.05)
			MvP <- dplyr::filter(MvP, interaction_qval > 0.05)
		}
		# Subclasses from tables, filter by WT results, best indication of significant fold change in response to treatment
		F_M <- dplyr::filter(FvM, qval.WT > 0.05)$target_id
		FM  <- dplyr::filter(FvM, qval.WT <= 0.05 & b<=0)$target_id
		MF  <- dplyr::filter(FvM, qval.WT <= 0.05 & b>0)$target_id
		F_P <- dplyr::filter(FvP, qval.WT > 0.05)$target_id
		FP  <- dplyr::filter(FvP, qval.WT <= 0.05 & b<=0)$target_id
		PF  <- dplyr::filter(FvP, qval.WT <= 0.05 & b>0)$target_id
		M_P <- dplyr::filter(MvP, qval.WT > 0.05)$target_id
		MP  <- dplyr::filter(MvP, qval.WT <= 0.05 & b<=0)$target_id
		PM  <- dplyr::filter(MvP, qval.WT <= 0.05 & b>0)$target_id
		# Intersection defines classes
		expression_class_list <- list(
			MFP = Reduce(intersect, list(MF, FP, MP)),
			PFM = Reduce(intersect, list(PF, FM, PM)),
			PM_F = Reduce(intersect, list(PM, PF, F_M)),
			F_MP = Reduce(intersect, list(F_M, MP, FP)),
			MF_P = Reduce(intersect, list(MF, MP, F_P)),
			F_PM = Reduce(intersect, list(FM, PM, F_P)),
			MPF = Reduce(intersect, list(MP, PF, MF)),
			PMF = Reduce(intersect, list(PM, MF, PF)),
			P_MF = Reduce(intersect, list(M_P, MF, PF)),
			FMP = Reduce(intersect, list(FM, MP, FP)),
			FPM = Reduce(intersect, list(FP, PM, FM)),
			FM_P = Reduce(intersect, list(FM, M_P, FP))
		)
		# Write to files
		output_path <- paste(
				results_path,
				g,
				'expression_classes',
				sep='/'
			)
		dir.create(output_path, recursive=TRUE)
		expression_classes_counts <- data.frame(lengths(expression_class_list))
		colnames(expression_classes_counts) <- 'counts'
 		expression_classes_counts <- cbind(expression_classes_counts, expression_classes)
		write.table(
			expression_classes_counts,
			file =  paste0(output_path, '/counts_per_class.txt'),
			row.names = FALSE,
			sep='\t'
		)
		print(expression_classes_counts)
		for (i in seq_along(expression_class_list)) {
			if (length(expression_class_list[i]) != 0) {
				write.table(
					expression_class_list[i],
					file = paste(
						output_path,
						paste0(
							expression_classes$code[i],
							'_',
							expression_classes$file_annotation[i],
							'.txt'
						),
						sep='/'
					),
					row.names=FALSE,
					col.names=FALSE,
					quote=FALSE,
					sep="\t"
				)
			}
		}
	}
}

