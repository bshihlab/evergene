#### Analyse button (process data) and associated dropdown menus
analyseButton <- function(id) {
	ns <- NS(id)
	tagList(
		textOutput(ns("analyse_message")),
		actionButton(ns("button_compute"), "Analyse", style='font-size:120%') 	
	)
}

analyseButtonServer <- function(id, input_cancer_type, input_gene_string, user_input_show, upload_count, upload_annotation, upload_gene, upload_survival, upload_log_transform, upload_raw_count, rv) {
  moduleServer(
	id,
	function(input, output, session) {
		observeEvent(input$button_compute, {
			validate(need(isTruthy(input_cancer_type()) || (isTruthy(upload_count()) & isTruthy(upload_annotation())), 
					message = "Please select at least 1 cancer project or upload data & annotation"))
			## Validate inputs -------------------------------------
			# Show warning messages if missing inputs
			
			#iv <- InputValidator$new()
			#iv$add_rule("gene_string", sv_required())
			#iv$add_rule("input_cancer_type",  sv_required())
			#iv$add_rule("gene_string", function(x){length(input_format(x)[["formatted_gene_list_named"]])>0})
			#iv$enable()
			# if user input upload menu is showing and there are uploaded annotation/count
			if(user_input_show() & (isTruthy(upload_count()) & isTruthy(upload_annotation())) ){
				processed_user_input <- userInput_load(in_count = upload_count(),
												in_annotation = upload_annotation(),
												in_gene = upload_gene(),
												in_survival = upload_survival(),
												input_gene_string = input_gene_string(),
												in_log_transform = upload_log_transform(),
												in_raw_count = upload_raw_count())
				
				if(isTruthy(processed_user_input[['pca']])){
					rv[["dataSource"]] = "User input"
					rv[["showPanel"]] <- TRUE
					rv[["exprs_unit"]] = "" 
					rv[["input_cancer_type"]] <- c("User input" = "User input")
					rv[["max_pc"]] <- processed_user_input[["max_pc"]]
					rv[["formatted_gene_list"]] <-  c(processed_user_input[["formatted_gene_list"]])
					rv[["pca"]] = list("User input" = processed_user_input[["pca"]])
					rv[["scree"]] = list("User input" = processed_user_input[["scree"]])
					rv[["survival"]] = list("User input" = processed_user_input[["survival"]]) 
					rv[["top_contributions"]] = list("User input" = processed_user_input[["top_contributions"]])
					rv[["dropdown_pca"]] = list("User input" = processed_user_input[["dropdown_pca"]])
					rv[["dropdown_survival"]] = list("User input" = processed_user_input[["dropdown_survival"]]) 
					rv[["dropdown_correlation"]] = list("User input" = processed_user_input[["dropdown_correlation"]])
					rv[["display_message"]] <- processed_user_input[["message"]]
					
					# if there are more than 2 genes, average to generate gene geneset
					if(length(rv[["formatted_gene_list"]] ) > 1 ){
						rv[["formatted_gene_list"]] <- c(rv[["formatted_gene_list"]], geneset_dropdown)
					}

				} else {
					rv[["dataSource"]] = "User input"
					rv[["showPanel"]] <- FALSE
					rv[["display_message"]] <- processed_user_input[["message"]]
				}
			} else {
				# Lable selected cancer types with readable names
				formatted_input_cancer_type <- input_cancer_type()
				names(formatted_input_cancer_type) <- list_of_cancer_types_rev[input_cancer_type()]
	
				# format gene names
				formatting_output <- input_format(input_gene_string(), input_cancer_type())
				formatted_gene_list <- formatting_output[["formatted_gene_list_named"]]
				
				# If none of the entered genes match the database, add TP63 as a default placeholder and show warning message.
				if(length(formatted_gene_list) == 0){ 
						formatted_gene_list <- c("Input gene not found: TP63 as placeholder" = "ENSG00000073282.14")
				}
				rv[["display_message"]] <- "Gene input not found in the dataset. Using TP63 as a placeholder."
					

				## Initiate the reactive values 
				rv[["dataSource"]] = "TCGA"
				rv[["showPanel"]] <- TRUE
				rv[["exprs_unit"]] = "TPM" 
				rv[["input_cancer_type"]] <- formatted_input_cancer_type
				rv[["formatted_gene_list"]] <- formatted_gene_list
				rv[["pca"]] = list()
				rv[["scree"]] = list()
				rv[["survival"]] = list()
				rv[["top_contributions"]] = list()
				rv[["dropdown_pca"]] = list()
				rv[["dropdown_survival"]] = list()
				rv[["dropdown_correlation"]] = list()
				rv[["display_message"]] = ""
				
				
				# if there are more than 2 genes, average to generate gene geneset
				if(length(formatted_gene_list) > 1 ){
					rv[["formatted_gene_list"]] <- c(rv[["formatted_gene_list"]], geneset_dropdown)
				}
				print(rv[["dropdown_gene"]])
				## Indicate computing is in progress
				# Use a value to indicate that compute is in progress and do not update reactive plots while this is in progress
				# This is to go around the issues that when users select a new set of cancers/gene while already having plot display, 
				# it would cause errors due to some values being midway through processing
				rv$computingInProgress <- TRUE

	

				## read in data
				# Wrap the read in process within a progress bar
				withProgress(message = 'Reading in data', value = 0, {
				  # Number of times we'll go through the loop
					progress_idx <- 1
					for(current_cancer_type in rv$input_cancer_type){
						# Get the full cancer name
						current_cancer_full_name <- list_of_cancer_types_rev[current_cancer_type]
						# Progress bar
						incProgress(1/length(rv$input_cancer_type), detail = paste("Generating part", progress_idx))
						progress_idx <- progress_idx+1
						
						# Read in all required data
						# top contributions
						rv[["top_contributions"]][[current_cancer_type]] <- read_fst(paste0("data/top_contributions/", current_cancer_type, ".fst"))
						

						# Only read in required genes for tpm. Tpm is currently organised so columns refer to genes
						# tpm required by both analysis
						tpm_data <- read_fst(paste0("data/tpm/", current_cancer_type, ".fst"), columns = c("barcode", formatted_gene_list))
						row.names(tpm_data) <- tpm_data$barcode

						#### Add in value for gene geneset when there is more than one gene
						if(length(formatted_gene_list) > 1 ){
							tpm_geneset_calculate <- log2(data.matrix(tpm_data[,formatted_gene_list] + 0.1))
							tpm_data$geneset <- 2^(rowMeans(tpm_geneset_calculate))
						}
				  
						# Data required for PCA
						# Annotate sample-PCA data with TPM
						sample_annotation <- read_fst(paste0("data/sample_annotation/", current_cancer_type, ".fst"))
						pca_data <- merge(sample_annotation, tpm_data, by="barcode")
						pca_data$PFI_status <- factor(pca_data$PFI, levels=c(0,1))
						levels(pca_data$PFI_status) <- c("Yes", "No")
						#pca_data$tumour_stage_shortened <- NULL # remove tumour staging
						rv[["pca"]][[current_cancer_type]] <- pca_data
						rv[["scree"]][[current_cancer_type]] <- read_fst(paste0("data/scree/", current_cancer_type, ".fst"))
						
						
						# Data required for survival analysis
						#tpm_data$barcode = NULL
						#rv[["tpm"]][[current_cancer_type]] <- tpm_data




																			
						## Organise plot data for survival
						# Remove duplicates from the same patient
						labelled_survival_data <- merge(sample_annotation, tpm_data, by.x="barcode", by.y=0)
						labelled_survival_data <- labelled_survival_data[!duplicated(labelled_survival_data$bcr_patient_barcode),]
						# Remove samples with no information on surv_status (OS or PFI) or surv_time (OS.time or PFI.time)
						# Store plot data for survival
						rv[["survival"]][[current_cancer_type]] <- labelled_survival_data

						## Organise the dropdown menu for each analysis
						dropdown_pca <- add_sampleAnnot_pca[add_sampleAnnot_pca %in% colnames(rv[["pca"]][[current_cancer_type]])]
						dropdown_corr <- c(rv[["formatted_gene_list"]], as.character(dropdown_pca[dropdown_pca %in% annotation_numeric]), pc_names)
						names(dropdown_corr) <- c(names(rv[["formatted_gene_list"]]), names(dropdown_pca)[dropdown_pca %in% annotation_numeric], names(pc_names))
						dropdown_survival <- add_sampleAnnot_survival[as.character(add_sampleAnnot_survival) %in% colnames(rv[["pca"]][[current_cancer_type]])]

						rv[["dropdown_pca"]][[current_cancer_type]]  <- dropdown_pca
						rv[["dropdown_correlation"]][[current_cancer_type]]  <- factor(dropdown_corr, levels=dropdown_corr)
						rv[["dropdown_survival"]][[current_cancer_type]]  <- dropdown_survival
					}
				})
			}

			
			# All compute has been finished and reactive plots can now be updated
			rv$computingInProgress <- FALSE
		})
		output$analyse_message <- renderText({
					validate(need((isTruthy(input_cancer_type()) || (isTruthy(upload_count()) & isTruthy(upload_annotation()))), message = "Please select at least 1 cancer project or upload data & annotation."),
							#need(isTruthy(input_gene_string()), message = "Please enter at least 1 gene.")
							)
					rv[["display_message"]]})
	})
}