#### This generates the pre-processed PCA analysis on all projects
# Date downloaded from TCGA using TCGAbiolink on 1-August-2023
#0 Set up----------------------------------------------------------------------

#0 Settings -----------------------------------------------------------------------
working_dir <- "D:/shared/nwcr_data/"
working_dir <- "C:/Users/shihb/OneDrive - Lancaster University/work/project/2023/20230816_nwcr/nwcr_data"
input_dir_tgca <- "raw"
input_fp_stemness <- "additional_annotation/stemness.csv"
input_fp_patientAnnotation <- "additional_annotation/patient_annotation.csv"
input_dir_subtype <- "additional_annotation/subtype.csv"
input_fp_event_recommendation <- "additional_annotation/TCGA-CDR_recommandation.csv"
output_dir_pca_scree <- "data/scree"
output_dir_pca_scree_contribQual <- "data/scree_contribQual"
output_dir_tpm <- "data/tpm"
output_dir_event_recommendation <- "data/event_recommendation"
output_dir_sampleAnnotation <- "data/sample_annotation"
output_dir_clinical <- "data/clinical"
patient_annotation_keep <- c("gender", "race", "vital_status", "ajcc_pathologic_tumor_stage", "clinical_stage", "histological_type", "age_at_initial_pathologic_diagnosis", "OS", "OS.time", "PFI", "PFI.time", "DSS", "DSS.time", "DFI", "DFI.time", "treatment_outcome_first_course", "new_tumor_event_dx_days_to", "new_tumor_event_type")
num_pcs <- 20

setwd(working_dir)

dir.create("data", showWarnings=FALSE)
dir.create(output_dir_tpm, showWarnings=FALSE)
dir.create(output_dir_sampleAnnotation, showWarnings=FALSE)
dir.create(output_dir_clinical, showWarnings=FALSE)
dir.create(output_dir_pca_scree, showWarnings=FALSE)
dir.create(output_dir_pca_scree_contribQual, showWarnings=FALSE)
dir.create(output_dir_event_recommendation, showWarnings=FALSE)

source("script/load_libraries.R")
source("script/load_variables.R")
library(edgeR)
library(SummarizedExperiment)
library(TCGAbiolinks)

stemness_df <- read.csv(input_fp_stemness)
patient_annotation <- read.csv(input_fp_patientAnnotation)
patient_annotation[patient_annotation == "#N/A"] <- NA
patient_annotation[patient_annotation == "[Unknown]"] <- NA
patient_annotation[patient_annotation == "[Not Available]"] <- NA
event_recommendation <- read.csv(input_fp_event_recommendation)

## read in immune subtype and tcga subtype information
immune_subtype <- read.csv("additional_annotation/immune_subtype.csv")
tcga_subtype_shortened <- read.csv("additional_annotation/tcga_subtype_shortened.csv")
# organise shortened subtype
immune_subtype <- merge(immune_subtype, tcga_subtype_shortened, by="tcga_subtype", all.x=TRUE)
immune_subtype$tcga_subtype_shortened <- ifelse(immune_subtype$tcga_subtype_shortened == "NA", NA, immune_subtype$tcga_subtype_shortened)
immune_subtype$tcga_subtype_shortened <- ifelse(immune_subtype$tcga_subtype_shortened == "", NA, immune_subtype$tcga_subtype_shortened)

# organise the event annotation
event_recommendation$display_text <- ifelse(event_recommendation$caution > 0, 
                                            paste(#event_recommendation$recommendation,
                                              "TCGA-CDR recommendation caution: ",    
                                              event_recommendation$caution_note), "" )
event_recommendation$cancer_type <- paste0("TCGA-", event_recommendation$cancer_type)
event_recommendation <- event_recommendation[,c("cancer_type", "event", "display_text")]
write_fst(event_recommendation, paste0(output_dir_event_recommendation, "/TCGA_CDR_recommandation.fst"))


# Shortened text for histological type
histological_type_shortened <- read.csv("additional_annotation/histological_type_shortened.csv")


extra_info <- c("figo_stage", "paper_DIFFERENTIATION", "ajcc_pathologic_stage", "race", "gender")
all_projects <- list.files(input_dir_tgca)
all_projects <- gsub(".rds", "", all_projects)

sample_annotation_list <- list()

#1 Loop through the cancer projects -------------------------------------------------------------------
for(current_cancer_type in list_of_cancer_types){
	print(current_cancer_type) # Prints the names of the files that are being accessed in the loop

	#2 current file paths -----------------------------------------------------------------------
	current_in_fp <- paste0(input_dir_tgca, "/", current_cancer_type, ".rds")
	current_out_fp_pca_scree <- paste0(output_dir_pca_scree, "/", current_cancer_type, ".fst")
	current_out_fp_pca_scree_contribQual <- paste0(output_dir_pca_scree_contribQual, "/", current_cancer_type, ".fst")
	current_out_fp_tpm <- paste0(output_dir_tpm, "/", current_cancer_type, ".fst")
	current_out_fp_sampleAnnotation <- paste0(output_dir_sampleAnnotation, "/", current_cancer_type, ".fst")
	current_out_fp_clinical <- paste0(output_dir_clinical, "/", current_cancer_type, ".fst")
	
	current_summarised_exp <- readRDS(current_in_fp)  


	#3 Organise count data for PCA -----------------------------------------------------------------------
	count_data <- assay(current_summarised_exp, "unstranded") # gene count matrix

	# Normalise data
	x <- DGEList(counts = count_data)
	keep.exprs <- filterByExpr(x, min.prop=0.2)	# genes expressed in at least 20% of the population
	x <- x[keep.exprs,, keep.lib.sizes=FALSE]
	x <- calcNormFactors(x, method = "TMM")

	lcpm_filtered <- cpm(x, log=TRUE)
	

	#4 PCA------------------------------------------------------------------------
	set.seed(1)
	count_pca <- prcomp(t(lcpm_filtered), scale=TRUE) # Carry out PCA on transformed matrix of count data.

	# Making scree plots----------
	# Foratting data: shorten, transform, merge
	pca_summary <- summary(count_pca)$importance[2,] * 100 # Showing the contribution of each position, to represent its importance.
	contribution_dataframe <- data.frame(PC = paste0("PC", 1:num_pcs), contribution = pca_summary[1:num_pcs])

	
	# Reorganise to calculate the strength of contribution for each gene.
	#### edit this
	# Results for Variables
	res.var <- get_pca_var(count_pca)
	gene_contribution <- res.var$contrib			# Gene contributions to the PCs
	gene_contribution_cos2 <- res.var$cos2			# Quality of representation 
	
	gene_contribution <- res.var$contrib[,1:num_pcs]
	# add the (+/-)sign for each gene to contribution
	gene_contribution <- as.data.frame(sign(res.var$coord[,1:num_pcs] )) * as.data.frame(gene_contribution)
	gene_contribution <- t(gene_contribution)
	gene_contribution_cos2 <- t(gene_contribution_cos2[,1:num_pcs])
	
	rownames(gene_contribution) <- paste0("PC", 1:num_pcs)
	rownames(gene_contribution_cos2) <- paste0("PC", 1:num_pcs)
	
	scree_df <- merge(contribution_dataframe, gene_contribution, by.x="PC", by.y=0)
	scree_df_cos2 <- merge(contribution_dataframe, gene_contribution_cos2, by.x="PC", by.y=0)
	write_fst(scree_df, current_out_fp_pca_scree, compress=100)
	write_fst(scree_df_cos2, current_out_fp_pca_scree_contribQual, compress=100)
	


	# The x-/y-coordinates to plot for top 10 PC
	pca_positions <- as.data.frame(count_pca$x)
	pca_positions <- pca_positions[,1:num_pcs]


	#5 Trim SummarizedExperiment data so it only keeps TPM and a subset of the sample annotation ----------------
	# Also add the PC coordinates for each sample for the first 10 PCs
	sample_annotation <- colData(current_summarised_exp)
	sample_annotation_trimmed <- data.frame(barcode = sample_annotation$barcode,
									sample = sample_annotation$sample, 
									sample_type = sample_annotation$sample_type)
									
	# add additional information if they exist
	sample_annotation_trimmed <- merge(sample_annotation_trimmed, stemness_df, by.x="barcode", by.y="TCGAlong.id", all.x=TRUE)
	#sample_annotation_trimmed <- merge(sample_annotation_trimmed, subtype_df[,c("barcode", "subtype")], by.x="barcode", by.y="barcode", all.x=TRUE) # using subtype from pancancer immune landscape paper instead
	sample_annotation_trimmed <- merge(sample_annotation_trimmed, pca_positions, by.x="barcode", by.y=0, all.x=TRUE)
	original_order <- sample_annotation_trimmed$barcode

	# merge with patient annotation
	sample_annotation_trimmed$bcr_patient_barcode <- substr(sample_annotation_trimmed$barcode, 1, 12)
	sample_annotation_trimmed <- merge(sample_annotation_trimmed, 
				patient_annotation[,c("bcr_patient_barcode", patient_annotation_keep)],
				by="bcr_patient_barcode", all.x=TRUE)
					
	row.names(sample_annotation_trimmed) <- sample_annotation_trimmed$barcode
	sample_annotation_trimmed <- sample_annotation_trimmed[original_order,c("barcode", colnames(sample_annotation_trimmed)[!colnames(sample_annotation_trimmed) %in% "barcode"])]
	sample_annotation_trimmed$race <- tolower(sample_annotation_trimmed$race)
	sample_annotation_trimmed$gender <- tolower(sample_annotation_trimmed$gender)
	sample_annotation_trimmed$clinical_stage <- gsub("\\[Not Applicable\\]", NA, sample_annotation_trimmed$clinical_stage )
	sample_annotation_trimmed$ajcc_pathologic_tumor_stage <- gsub("\\[Not Applicable\\]", NA, sample_annotation_trimmed$ajcc_pathologic_tumor_stage )
	sample_annotation_trimmed$tumour_stage <- ifelse(is.na(sample_annotation_trimmed$ajcc_pathologic_tumor_stage), sample_annotation_trimmed$clinical_stage, sample_annotation_trimmed$ajcc_pathologic_tumor_stage) 
	sample_annotation_trimmed$ajcc_pathologic_tumor_stage <- gsub("\\[\\]", "", sample_annotation_trimmed$ajcc_pathologic_tumor_stage )
	sample_annotation_trimmed$clinical_stage <- NULL
	sample_annotation_trimmed$ajcc_pathologic_tumor_stage <- NULL

	# merge to label immune subtype 
	sample_annotation_trimmed <- merge(sample_annotation_trimmed, immune_subtype[,c("bcr_patient_barcode", "immune_subtype", "tcga_subtype_shortened")], by="bcr_patient_barcode", all.x=TRUE)
	
	# remove clinical stage, tumour, stage, histological_type, or ethnicity if all samples are identical
	for(annotation in colnames(sample_annotation_trimmed)){
		current_unique_vals <- unique(sample_annotation_trimmed[[annotation]])
		current_unique_vals_no_na <- current_unique_vals[!is.na(current_unique_vals)]
		if(length(current_unique_vals_no_na) <= 1 ){
			sample_annotation_trimmed <- sample_annotation_trimmed[, colnames(sample_annotation_trimmed)[!colnames(sample_annotation_trimmed) %in% annotation]] 
		}	
	}

	# Shorten tumour stage and turn it into numerics
	if("tumour_stage" %in% colnames(sample_annotation_trimmed)){
		sample_annotation_trimmed$tumour_stage_shortened <- ifelse(sample_annotation_trimmed$tumour_stage %in% c("[Discrepancy]", "[Unknown]", "Stage X", "I/II NOS"), NA, sample_annotation_trimmed$tumour_stage)
		sample_annotation_trimmed$tumour_stage_shortened <- tolower(sample_annotation_trimmed$tumour_stage_shortened )
		sample_annotation_trimmed$tumour_stage_shortened <- gsub("stage|a|b|c|s| |1|2", "", sample_annotation_trimmed$tumour_stage_shortened)
		sample_annotation_trimmed$tumour_stage_shortened <- as.numeric(as.roman(sample_annotation_trimmed$tumour_stage_shortened))
	}
	# Shorten histological_type
	if("histological_type" %in% colnames(sample_annotation_trimmed)){
	  sample_annotation_trimmed$histological_type_shortened = NULL
		sample_annotation_trimmed <- merge(sample_annotation_trimmed, histological_type_shortened[,c("histological_type", "histological_type_shortened")], by="histological_type", all.x=TRUE)
		sample_annotation_trimmed$histological_type_shortened <- ifelse(sample_annotation_trimmed$histological_type_shortened == "NA", NA, sample_annotation_trimmed$histological_type_shortened )
		sample_annotation_trimmed$histological_type_shortened <- ifelse(sample_annotation_trimmed$histological_type_shortened == "", NA, sample_annotation_trimmed$histological_type_shortened )
		sample_annotation_trimmed$histological_type_shortened <- ifelse(sample_annotation_trimmed$histological_type_shortened == "[Discrepancy]", NA, sample_annotation_trimmed$histological_type_shortened )
		sample_annotation_trimmed$histological_type_shortened <- gsub('(.{1,20})(\\s|$)', '\\1\n', sample_annotation_trimmed$histological_type_shortened)
		sample_annotation_trimmed$histological_type_shortened <- gsub('\\\n$', "", sample_annotation_trimmed$histological_type_shortened)
	}
	# Shorten race into shorter groups names
	if("race" %in% colnames(sample_annotation_trimmed)){
		sample_annotation_trimmed$race <- ifelse(sample_annotation_trimmed$race %in% c("[not evaluated]", "[Unknown]"), NA, sample_annotation_trimmed$race)
		sample_annotation_trimmed$race_shortened <- sample_annotation_trimmed$race
		sample_annotation_trimmed$race_shortened <- ifelse(sample_annotation_trimmed$race_shortened == "black or african american", "black", sample_annotation_trimmed$race_shortened)
		sample_annotation_trimmed$race_shortened <- ifelse(sample_annotation_trimmed$race_shortened == "american indian or alaska native", "other", sample_annotation_trimmed$race_shortened)
		sample_annotation_trimmed$race_shortened <- ifelse(sample_annotation_trimmed$race_shortened == "native hawaiian or other pacific islander", "other", sample_annotation_trimmed$race_shortened)
	}
	if("new_tumor_event_type" %in% colnames(sample_annotation_trimmed)){
		sample_annotation_trimmed$new_tumor_event_type_shortened <- sample_annotation_trimmed$new_tumor_event_type
		sample_annotation_trimmed$new_tumor_event_type_shortened <- ifelse(sample_annotation_trimmed$new_tumor_event_type_shortened == "Distant Metastasis|Distant Metastasis|Regional lymph node", "Regional lymph node|Distant Metastasis", sample_annotation_trimmed$new_tumor_event_type_shortened)
		sample_annotation_trimmed$new_tumor_event_type_shortened <- ifelse(sample_annotation_trimmed$new_tumor_event_type_shortened == "Distant Metastasis|Regional lymph node", "Regional lymph node|Distant Metastasis", sample_annotation_trimmed$new_tumor_event_type_shortened)
		sample_annotation_trimmed$new_tumor_event_type_shortened <- ifelse(sample_annotation_trimmed$new_tumor_event_type_shortened == "Regional lymph node|Distant Metastasis|Distant Metastasis", "Regional lymph node|Distant Metastasis", sample_annotation_trimmed$new_tumor_event_type_shortened)
	}
	if("treatment_outcome_first_course" %in% colnames(sample_annotation_trimmed)){
		sample_annotation_trimmed$treatment_outcome_first_course_shortened <- sample_annotation_trimmed$treatment_outcome_first_course
		sample_annotation_trimmed$treatment_outcome_first_course_shortened <- ifelse(sample_annotation_trimmed$treatment_outcome_first_course_shortened == "Normalization of Tumor Markers, but Residual Tumor Mass", "Other", sample_annotation_trimmed$treatment_outcome_first_course_shortened)
		sample_annotation_trimmed$treatment_outcome_first_course_shortened <- ifelse(sample_annotation_trimmed$treatment_outcome_first_course_shortened == "[Discrepancy]", NA, sample_annotation_trimmed$treatment_outcome_first_course_shortened)
		sample_annotation_trimmed$treatment_outcome_first_course_shortened <- ifelse(sample_annotation_trimmed$treatment_outcome_first_course_shortened == "[Not Applicable]", NA, sample_annotation_trimmed$treatment_outcome_first_course_shortened)
		sample_annotation_trimmed$treatment_outcome_first_course_shortened <- ifelse(sample_annotation_trimmed$treatment_outcome_first_course_shortened == "[Not Evaluated]", NA, sample_annotation_trimmed$treatment_outcome_first_course_shortened)
		sample_annotation_trimmed$treatment_outcome_first_course_shortened <- ifelse(sample_annotation_trimmed$treatment_outcome_first_course_shortened == "[Unknown]", NA, sample_annotation_trimmed$treatment_outcome_first_course_shortened)
		sample_annotation_trimmed$treatment_outcome_first_course_shortened <- ifelse(sample_annotation_trimmed$treatment_outcome_first_course_shortened == "No Measureable Tumor or Tumor Markers", "No Measureable Tumor", sample_annotation_trimmed$treatment_outcome_first_course_shortened)
	}
# 	no longer required as query has limited to primary tumour
#	add annotation on whether the sample is tumour or normal
#	NT_barcode <- TCGAquery_SampleTypes(sample_annotation_trimmed$barcode, typesample = c("NT"))
#	TP_barcode <- TCGAquery_SampleTypes(sample_annotation_trimmed$barcode, typesample = c("TP"))
#	sample_annotation_trimmed$status <- ifelse(sample_annotation_trimmed$barcode %in% NT_barcode, "normal", NA)
#	sample_annotation_trimmed$status <- ifelse(sample_annotation_trimmed$barcode %in% TP_barcode, "tumour", sample_annotation_trimmed$barcode %in% TP_barcode)
	
	# rename some sample types
#	sample_annotation_trimmed$sample_type <- gsub("Primary Blood Derived Cancer - Peripheral Blood", "Primary Blood \nDerived Cancer", sample_annotation_trimmed$sample_type)
#	sample_annotation_trimmed$sample_type <- gsub("Additional - ", "", sample_annotation_trimmed$sample_type)
#	sample_annotation_trimmed$sample_type <- gsub("Additional - ", "", sample_annotation_trimmed$sample_type)
	
	write_fst(sample_annotation_trimmed, current_out_fp_sampleAnnotation)
	

	# Save the tpm as dataframe. Change it so gene_id is column name
	tpm_out <- as.data.frame(t(assay(current_summarised_exp, "tpm_unstrand")))
	tpm_out$barcode <- row.names(tpm_out)
	write_fst(tpm_out, current_out_fp_tpm, compress=100)

	#6 Clinical data 
	data_clinic <- GDCquery_clinic(current_cancer_type, "clinical")
	write_fst(data_clinic, current_out_fp_clinical)
	
} 
 



