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

#Additional sample annotations
add_sampleAnnot_pca <- c(
  "Subtypes: molecular subtype" =  "tcga_subtype_shortened", 
  "Subtypes: immune type" = "immune_subtype",
  "Subtypes: histological type" = "histological_type_shortened",
#  "Tumour: stage" =  "tumour_stage_shortened", 
  "Tumour: treatment outcome (first course)" = "treatment_outcome_first_course_shortened",
  "Tumour: new tumour event type" = "new_tumor_event_type_shortened",
  "Molecular: stemness (mRNAsi)" = "mRNAsi",
  "Molecular: stemness(EREG.mRNAsi)" = "EREG.mRNAsi",
  "Demographic: age" = "age_at_initial_pathologic_diagnosis",
  "Demographic: race" = "race_shortened",
  "Demographic: gender" =  "gender", 
  "Status: vital status" = "vital_status",
  "Time: days till death" = "OS.time_dead",
  "Time: days till new tumour event" = "new_tumor_event_dx_days_to"
)

add_sampleAnnot_survival <- c(
  "Status: vital status" = "vital_status",
  "Time: days tracked" = "OS.time",
  "Time: days till death" = "OS.time_dead",
  "Time: days till new tumour event" = "new_tumor_event_dx_days_to",
  "Subtypes: molecular subtype" =  "tcga_subtype_shortened", 
  "Subtypes: immune type" = "immune_subtype",
  "Subtypes: histological type" = "histological_type_shortened",
#  "Tumour: stage" =  "tumour_stage_shortened",
  "Tumour: treatment outcome (first course)" = "treatment_outcome_first_course_shortened",
  "Tumour: new tumour event type" = "new_tumor_event_type_shortened",
  "Molecular: stemness (mRNAsi)" = "mRNAsi",
  "Molecular: stemness(EREG.mRNAsi)" = "EREG.mRNAsi",
  "Demographic: age" = "age_at_initial_pathologic_diagnosis",
  "Demographic: race" = "race_shortened",
  "Demographic: gender" =  "gender"
)

# tumour: staging is ajcc_pathologic_tumor_stage or clinical_stage
# list numeric annotations
annotation_numeric <- c("OS.time_dead", "OS.time", "new_tumor_event_dx_days_to", "age_at_initial_pathologic_diagnosis", 
	"tumour_stage_shortened", "mRNAsi", "EREG.mRNAsi")
annotation_jitter <- c("race_shortened", "gender", "immune_subtype", "vital_status", 
	"tumour_stage_shortened", "histological_type_shortened", "tcga_subtype_shortened", 
	"new_tumor_event_type_shortened", "treatment_outcome_first_course_shortened")


#### Compute
out_dir_pca <- "evergene_pca"
out_dir_survival <- "evergene_survival"
maximum_gene_project_count <- 500
top_gene_contribution <- 20
analysis_computed <- FALSE


#### Help messages
help_gene_input <- c("Genes formated as gene name or Ensembl ID, separately by spaces, commas or new lines.", 
	"Input genes not identified will be indicated at the bottom of the side panel")
help_cancer_input <- c("Tick individual cancers of interest from the dropdown menu. You can click on individual cancer types to select one or more cancer types, or the <b>select all</b> or <b>unselect all</b> buttons to select/unselect all cancer types. 
					Selected cancer types have a tick next to it. To exit the dropdown meanu, click on anywhere outside the menu.", 
					paste0('Primary tumours from a total of 32 cancer types are available. RNA sequencing are derived from the ',
						as.character(tags$a("TCGA Harmonized data", href= "https://docs.gdc.cancer.gov/Encyclopedia/pages/Harmonized_Data/", target="_blank")),
						' through ', 
						as.character(tags$a("TCGAbiolinks", href= "doi:10.1093/nar/gkv1507", target="_blank")), 
						". Details on the RNA processing pipelines can be found in ", 
						as.character(tags$a("GDC documentation, mRNA Analysis Pipeline ", href= "https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/", target="_blank")),
						". Clinical patient data are desrived from ", as.character(tags$a("PanCancer Atlas", href= "doi:10.1016/j.cell.2018.02.052", target="_blank")), "."))
help_example1 <- c("Description on Example 1") 
help_example2 <- c("Description on Example 2") 
help_example3 <- c("Description on Example 3") 

help_rowMenu <- c("<b>Gene</b>: dropdown menu showing the input genes used in the analysis.", 
					"<b>Cancer type</b>: dropdown menu showing the input cancer types used in the analysis.", 
					paste0("<b>PC (x-axis)</b>: the principal component (PC) plotted in the x-axis for PCA (principal component analysis) plot A and B. 
							PCs are generated from principal component analysis (PCA), which is a dimension reduction technique that is used 
							to simplify the complexity of data while retaining trends and patterns. You can find more information on PCA ",
							as.character(tags$a("here", href= "https://www.nature.com/articles/nmeth.4346", target="_blank")), "."),
					paste0("<b>PC (y-axis)</b>: the principal component plotted in the y-axis for PCA (principal component analysis) or survival analysis plot B. 
							PCs are generated from principal component analysis (PCA), which is a dimension reduction technique that is used 
							to simplify the complexity of data while retaining trends and patterns. You can find more information on PCA ",
							as.character(tags$a("here", href= "https://www.nature.com/articles/nmeth.4346", target="_blank")), "."))

help_gene_contribution <- c("") 
help_pca_results <- c("") 
help_survival_results <- c("") 

				

#### Colours
discrete_cols_exprs <- c("Low.exprs" = "#89CFF0", "Mid.exprs" = "#a3a6a7", "High.exprs" = "#cc210a")
discrete_cols_surv <- c(discrete_cols_exprs[1], discrete_cols_exprs[3])
discrete_cols <- c("#FF7276", "#00008B", "#33a02c", "#cab2d6", "#6a3d9a", "#b15928", "#ff7f00", "#a6cee3", "#b2df8a")




