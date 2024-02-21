### Set up ###
library(fts)
library(TCGAbiolinks)
library(SummarizedExperiment)


out_dir <- "nwcr_data/raw/"

### Providing a list of the cancer types
sample_count_threshold <- 80

cancer_codes <- c("TCGA-BRCA", "TCGA-THCA", "TCGA-UCEC", "TCGA-DLBC", "TCGA-COAD", 
				"TCGA-CESC", "TCGA-BLCA", "TCGA-CHOL", "TCGA-ESCA", "TCGA-ACC", 
				"TCGA-KICH", "TCGA-HNSC", "TCGA-LIHC", "TCGA-MESO", "TCGA-LAML", 
				"TCGA-KIRP", "TCGA-KIRC", "TCGA-GBM", "TCGA-LGG", "TCGA-SARC", 
				"TCGA-PCPG", "TCGA-READ", "TCGA-PAAD", "TCGA-LUAD", "TCGA-PRAD", 
				"TCGA-OV", "TCGA-LUSC", "TCGA-TGCT", "TCGA-THYM", "TCGA-UVM", 
				"TCGA-SKCM", "TCGA-UCS", "TCGA-STAD")

small_project <- list()

### Query and removing any <100 ###
cancer_codes_passed <- list()
for(x in cancer_codes){
	project <- x
	query_RNASeq <- GDCquery(project = project , 
						   data.category = "Transcriptome Profiling", 
						   data.type = "Gene Expression Quantification", 
						   workflow.type = "STAR - Counts", 
						   experimental.strategy = "RNA-Seq",
						   sample.type = c("Primary Tumor"))
	RNAseq_samples <- getResults(query_RNASeq)

	if(nrow(RNAseq_samples) > sample_count_threshold ){
	cancer_codes_passed[[x]] <- query_RNASeq
	print(nrow(RNAseq_samples))

	} else {
	small_project[[project]] <- nrow(RNAseq_samples)
	}

}


### Download ###
for (x in names(cancer_codes_passed)){
	GDCdownload(cancer_codes_passed[[x]]) 
	out <- GDCprepare(cancer_codes_passed[[x]])
	save_fst(out, paste0(out_dir, x, ".rds"))

}


