#### Background variables
list_of_cancer_types <- c(
#                          "Adrenocortical Carcinoma" = "TCGA-ACC",                                                                              
                          "Bladder Urothelial Carcinoma" = "TCGA-BLCA",                                                                              
                          "Breast Invasive Carcinoma" = "TCGA-BRCA",                                                                             
                          "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma" = "TCGA-CESC",                                      
#                          "Cholangiocarcinoma" = "TCGA-CHOL",                                                                                    
                          "Colon Adenocarcinoma" = "TCGA-COAD",                                                                                  
#                          "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma" = "TCGA-DLBC",                                                      
                          "Esophageal Carcinoma" = "TCGA-ESCA",                                                                                  
                          "Glioblastoma Multiforme" = "TCGA-GBM",                                                                               
                          "Head and Neck Squamous Cell Carcinoma" = "TCGA-HNSC",                                                
#                          "Kidney Chromophobe" = "TCGA-KICH",                                                                                    
                          "Kidney renal clear cell carcinoma" = "TCGA-KIRC",                                                                                    
                          "Kidney renal papillary cell carcinoma" = "TCGA-KIRP",                                                                                    
                          "Brain Lower Grade Glioma" = "TCGA-LGG",                                                                              
                          "Liver Hepatocellular Carcinoma" = "TCGA-LIHC",      
                          "Lung Adenocarcinoma" = "TCGA-LUAD",                                                                                   
                          "Lung Squamous Cell Carcinoma" = "TCGA-LUSC",                                                                        
                          "Mesothelioma" = "TCGA-MESO",                                                                     
                          "Ovarian Serous Cystadenocarcinoma" = "TCGA-OV",                                                                       
                          "Pheochromocytoma and Paraganglioma" = "TCGA-PAAD",                                                                    
#                          "Pheochromocytoma and Paraganglioma" = "TCGA-PCPG",  # removed due to high survival                                                                   
                          "Prostate Adenocarcinoma" = "TCGA-PRAD",                                                                               
                          "Rectum Adenocarcinoma" = "TCGA-READ",                                                                                 
                          "Sarcoma" = "TCGA-SARC",                                                                                               
                          "Skin Cutaneous Melanoma" = "TCGA-SKCM",                                                                               
                          "Stomach Adenocarcinoma" = "TCGA-STAD", 
                          "Testicular Germ Cell Tumors" = "TCGA-TGCT",                                                                           
                          "Thyroid Carcinoma" = "TCGA-THCA",                                                                                     
                          "Thymoma" = "TCGA-THYM",                                                                                               
                          "Uterine Corpus Endometrial Carcinoma" = "TCGA-UCEC",                                                                                               
#                          "Uterine Carcinosarcoma" = "TCGA-UCS",                                                                                
                          "Uveal Melanoma" = "TCGA-UVM")
list_of_cancer_types <- list_of_cancer_types[order(names(list_of_cancer_types))]
# Create a reversed name:value vector for list of cancer type
list_of_cancer_types_rev <- vector()
for(current_name in names(list_of_cancer_types)){
  current_code <- list_of_cancer_types[current_name]
  list_of_cancer_types_rev[current_code] = current_name
}

# list of Genes
geneset_dropdown <- c("Gene set" = "geneset")

#Additional sample annotations
add_sampleAnnot_pca <- c(
  "Subtypes: molecular subtype" =  "tcga_subtype_shortened", 
  "Subtypes: immune type" = "immune_subtype",
  "Subtypes: histological type" = "histological_type_shortened",
  "Tumour: stage" =  "tumour_stage_shortened", 
  "Tumour: treatment outcome (first course)" = "treatment_outcome_first_course_shortened",
  "Tumour: new tumour event type" = "new_tumor_event_type_shortened",
  "Molecular: stemness (mRNAsi)" = "mRNAsi",
  "Molecular: stemness (EREG.mRNAsi)" = "EREG.mRNAsi",
  "Demographic: age" = "age_at_initial_pathologic_diagnosis",
  "Demographic: race" = "race_shortened",
  "Demographic: gender" =  "gender", 
  "Event: vital status" = "vital_status",
  "Event: progression-free status"= "PFI_status",
  "Time: days till death" = "OS.time_dead",
  "Time: progression-free interval"= "PFI.time",
  "Time: days till new tumour event" = "new_tumor_event_dx_days_to"
)
add_sampleAnnot_pca <- factor(add_sampleAnnot_pca, levels=add_sampleAnnot_pca)

add_sampleAnnot_survival <- c(
  "Event: progression-free"= "PFI_status",
  "Event: disease-free" = "DFI_status", 
  "Event: disease-specific survival" = "DSS_status", 
  "Event: overall survival" = "vital_status", 
  "Time: progression-free interval"= "PFI.time",  
  "Time: days tracked" = "OS.time",
  "Time: days till death" = "OS.time_dead",
  "Time: days till new tumour event" = "new_tumor_event_dx_days_to",
  "Subtypes: molecular subtype" =  "tcga_subtype_shortened", 
  "Subtypes: immune type" = "immune_subtype",
  "Subtypes: histological type" = "histological_type_shortened",
  "Tumour: stage" =  "tumour_stage_shortened",
  "Tumour: treatment outcome (first course)" = "treatment_outcome_first_course_shortened",
  "Tumour: new tumour event type" = "new_tumor_event_type_shortened",
  "Molecular: stemness (mRNAsi)" = "mRNAsi",
  "Molecular: stemness (EREG.mRNAsi)" = "EREG.mRNAsi",
  "Demographic: age" = "age_at_initial_pathologic_diagnosis",
  "Demographic: race" = "race_shortened",
  "Demographic: gender" =  "gender"
)
add_sampleAnnot_survival <- factor(add_sampleAnnot_survival, levels=add_sampleAnnot_survival)

clinical_events_all <- c("Progression-free interval" = "PFI", "Overall survival"= "OS", "Disease-specific survivial"= "DSS", "Disease-free interval"= "DFI", "User input" = "user_input_event")

# tumour: staging is ajcc_pathologic_tumor_stage or clinical_stage
# list numeric annotations
annotation_numeric <- c("OS.time", "new_tumor_event_dx_days_to", "age_at_initial_pathologic_diagnosis", 
	"tumour_stage_shortened", "mRNAsi", "EREG.mRNAsi", "PFI.time", "DFI.time", "DSS.time")
annotation_jitter <- c("race_shortened", "gender", "immune_subtype", "vital_status", "PFI_status", "DSS_status", "DFI_status",
	"tumour_stage_shortened", "histological_type_shortened", "tcga_subtype_shortened", "treatment_outcome_first_course",
	"new_tumor_event_type_shortened", "treatment_outcome_first_course_shortened")

#### Read in TCGA gene annotation
tcga_gene_annotation <- read_fst("data/gene_annotation/gene_annotation.fst")

#### Compute
out_dir_pca <- "evergene_pca"
out_dir_survival <- "evergene_survival"
out_dir_correlation <- "evergene_correlation"
top_gene_contribution_num <- 20
num_pcs <- 20
analysis_computed <- FALSE
pc_names <- paste0("PC", 1:num_pcs)
names(pc_names) <- pc_names

#### Help messages
help_gene_input <- c(paste0("Input genes are used as annotation in PCA (performed on all genes) plots and for grouping in survival analysis (SA). Genes can be formatted as <b>gene names</b> or <b>Ensembl ID</b>, separately by <b>spaces</b>, <b>commas</b> or <b>new lines</b>. ", 
	"Input genes not identified will be indicated below the ANALYSE button."), 
	paste0("The maximum number of input is <b>100 gene-cancer combinations</b>. Numbers beyond that will automatically be filtered by limiting the input genes in the input order. ",
	"Large numbers of genes can take a few minutes for the loading and downloading processes."))
help_cancer_input <- c("Tick individual cancer types of interest from the dropdown menu (selected cancer types have a tick next to them). You can click on <b>multiple</b> cancer types individually or use the select all or unselect all buttons, or the <b>select all</b> or <b>unselect all</b> buttons to select/unselect all cancer types. 
					Selected cancer types have a tick next to it. To exit the dropdown menu, click on anywhere outside the menu.", 
					paste0('Primary tumours from a total of 26 cancer types are available. RNA sequencing data are derived from the ',
						as.character(tags$a("TCGA Harmonized data", href= "https://docs.gdc.cancer.gov/Encyclopedia/pages/Harmonized_Data/", target="_blank")),
						' through ', 
						as.character(tags$a("TCGAbiolinks", href= "doi:10.1093/nar/gkv1507", target="_blank")), 
						". Details on the RNA sequencing data processing pipelines can be found in ", 
						as.character(tags$a("GDC documentation, mRNA Analysis Pipeline ", href= "https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/", target="_blank")),
						". Clinical patient data are derived from ", as.character(tags$a("PanCancer Atlas", href= "doi:10.1016/j.cell.2018.02.052", target="_blank")), '".")',
						"The maximum number of input is <b>100 gene-cancer combinations</b>. Numbers beyond that will automatically be filtered by limiting the input genes in the input order."))
help_example1 <- c(paste0("TP63, SOX2, GATA4 and GATA6 all show strong contributions to Principle Component (PC) 1 , the PC that explains a large portion of the variations seen in the data (see Scree plot under Additional PCA plots). ", 
                          "The separation of the samples according to PC1 corresponds to its molecular subtype and histological subtype. ",  
                          as.character(tags$a("The Cancer Genome Atlas Research Network", href= "https://doi.org/10.1038/nature20805", target="_blank")), 
                          " has found TP63 and SOX2 to show frequent amplifications in squamous cell carcinoma, and GATA4 and GATA6 to be commonly amplified in adenocarcinomas.",
                          " In addition, the ethnic differences between the two subtypes can also been seen when changing the sample annotation to Demographic: race."),
                   paste0("PHB is the top gene contributing to PC2. PHB has been suggested to promote cancer cell proliferation in ESCC, gallbladder cancer, and bladder cancer (", 
                          as.character(tags$a("Young et al. 2018 ", href= "https://www.nature.com/articles/s41419-018-0661-3", target="_blank")),
                          ")."))
help_example2 <- c(paste0("Breast invasive carcinoma and key markers for its subtypes."))
help_example3 <- c("Description on Example 3") 

help_rowMenu <- c(	"<b>Cancer type</b>: drop-down menu showing the input cancer types used in the analysis.", 
					"<b>Gene/PC</b>: drop-down menu showing the input genes and PCs generated by the PCA analysis. gene set is the averaged value across all entered genes; the average is calculated as log average if it's TCGA data, custom input with count data or 'add log transformation' option.", 
					paste0("<b>Principal Component Axis</b>: Principal components (PCs) are generated from principal component analysis (PCA), ", 
							"which is a dimension reduction technique that is used to simplify the complexity of data while retaining trends and patterns. ",
							"You can find more information on PCA ",
							as.character(tags$a("here", href= "https://www.nature.com/articles/nmeth.4346", target="_blank")), ".)"),
					"<b>PC (x-axis)</b>: the PC plotted in the x-axis for PCA plots A and B.",
					"<b>PC (y-axis)</b>: the PC plotted in the y-axis for PCA plots A and B or survival analysis tab plot B.")

help_gene_contribution <- c(paste0("The % contribution for the selected gene in each principal component (PC).",  
						"Within each PC, the sum of the % contribution for all genes add up to 100%. The darker the colour, the stronger the contribution. ",
						"<b>Grey</b> indicates that the gene has low expression across samples for the given cancer type and was not included in the PCA."),
						"This graph is for helping users identify PCs they might be interested in.",
						paste0("Note that the earlier PCs explain larger portions of the variation within the data. Please see the <b>Scree plot</b> under 
						<b>Additional PCA plots</b> in the PCA tab for details on how much % variation each PC explains).") )
help_gene_contribution_heatmap <- c(paste0("Each column indicates a principal component (PC). The figure shows the top 20 genes with the highest % contribution in each PC. ",
									"Due to limitation on space, the users may need to hover over individual genes to read the gene name. ",
									"By hovering over the values, the users can also find out the % of contribution and whether the contribution is in the positive or negative direction. ",
									"Users may wish to compare the % contribution for their gene of interest to the % contribution for the top genes in the same PCs."),
									"Note that the earlier PC explain a higher proportion of the total variation within the data. Please refer to the Scree plot (plot K under Additional PCA plots) for the % variance explained by each PC.")


#### Colours
discrete_cols_exprs <- c("Low.exprs" = "#89CFF0", "Mid.exprs" = "#a3a6a7", "High.exprs" = "#cc210a")
discrete_cols_surv <- c(discrete_cols_exprs[1], discrete_cols_exprs[3])
discrete_cols <- c("#FF7276", "#00008B", "#33a02c", "#cab2d6", "#6a3d9a", "#b15928", "#ff7f00", "#a6cee3", "#b2df8a")


#### plot combination restriction
pca_download_combination_restriction <- 100
pca_plot_combination_restriction <- 100
maximum_gene_project_count <- 100


#### Survival plot annotation
clinical_event <- c("PFI" = "progression-free interval",
								"DFI" = "disease-free interval",
								"DSS" = "disease-specific survival",
								"OS" = "overall survival")
event_recommendation_all <- c("Progression-free interval" = "PFI", 
								"Disease-free interval" = "DFI",
								"Disease-specific survival" = "DSS",
								"Overall survival"= "OS")


organised_data <- list("tpm"=list(), 
					   "pca" = list(), 
					   "scree" = list(),
					   "survival" = list(),
					   "top_contributions" = list()) 