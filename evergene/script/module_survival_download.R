#### PCA: Download PDF
survivalDownloadPdf <- function(id, display_survival_event, survival_percent_lower, survival_percent_upper,  rv) {
  moduleServer(id,
	function(input, output, session) {
		ns <- session$ns
		observeEvent(
			input$action,
			{ showModal(
				modalDialog(
				  title = NULL,
				  h3("Download the file?", style = "text-align: center;"),
				  footer = tagList(
					downloadButton(
					  outputId = ns("download"),
					  label = "Yes"
					),
					modalButton("Cancel")
				  ),
				  easyClose = TRUE,
				  size = "m"
			))}
		)  # observeEvent
		output$download <- downloadHandler(
		filename = function() {"evergene_survival.pdf"}, 
		content = function(file) {
		on.exit(removeModal())
		shiny::withProgress(
			message = "Saving - this could take some time.",
			value = 0,
			{
				# Work out the genes to plot - only plot the first pca_download_combination_restriction gene-cancer combinations
				input_cancer_type <- rv[["input_cancer_type"]]
				max_gene_keep <-floor(pca_download_combination_restriction / length(input_cancer_type))
				if(length(rv[["formatted_gene_list"]]) > max_gene_keep){
				  plot_genes <- rv[["formatted_gene_list"]][1:max_gene_keep]
				  names(plot_genes) <- names(rv[["formatted_gene_list"]])[1:max_gene_keep]
				} else {
				  plot_genes <-rv[["formatted_gene_list"]] 
				}
				
				
				# open pdf file
				pdf(file, onefile = TRUE, width = 5, height = 5)
				progress_idx <- 1
				for(idx in 1:length(input_cancer_type)){
					current_cancer_type <- input_cancer_type[idx]
					current_cancer_full_name <- names(input_cancer_type)[idx]
					# Make plot for each gene
					for(gene in plot_genes){
						current_plot <- make_survival_plot(
							cancer_full_name = names(input_cancer_type)[input_cancer_type == current_cancer_type], 
							gene = gene,
							gene_full_name = names(plot_genes)[plot_genes == gene],
							labelled_data = rv[["survival"]][[current_cancer_type]],
							threshTop = survival_percent_upper(),
							threshBottom = survival_percent_lower(),
							event = display_survival_event())
						# print the survival plots
						print(current_plot[["surv_plot"]])
						# progress bar
						shiny::incProgress(progress_idx/(length(input_cancer_type) *length(plot_genes) ))
						progress_idx <- progress_idx + 1
					}
				}
				dev.off()
			}
			) # withProgress
			} # Content 
		) # downloadHandler
	} # function
  ) # moduleServer
} #function


		
#### PCA: Download zip
survivalDownloadZip <- function(id, display_survival_event, survival_percent_lower, survival_percent_upper, rv) {
  moduleServer(id,
	function(input, output, session) {
		ns <- session$ns
		observeEvent(
			input$action,
			{ showModal(
				modalDialog(
				  title = NULL,
				  h3("Download the file?", style = "text-align: center;"),
				  footer = tagList(
					downloadButton(
					  outputId = ns("download"),
					  label = "Yes"
					),
					modalButton("Cancel")
				  ),
				  easyClose = TRUE,
				  size = "m"
			))}
		)  # observeEvent
		output$download <- downloadHandler(
		filename = function(){"evergene_survival.zip"},
		content = function(file) {
			on.exit(removeModal())
			shiny::withProgress(
			message = "Saving - this could take some time.",
			value = 0,
			{
				# Work out the genes to plot - only plot the first pca_download_combination_restriction gene-cancer combinations
				input_cancer_type <- rv[["input_cancer_type"]]
				max_gene_keep <-floor(pca_download_combination_restriction / length(input_cancer_type))
				if(length(rv[["formatted_gene_list"]]) > max_gene_keep){
				  plot_genes <- rv[["formatted_gene_list"]][1:max_gene_keep]
				  names(plot_genes) <- names(rv[["formatted_gene_list"]])[1:max_gene_keep]
				} else {
				  plot_genes <-rv[["formatted_gene_list"]] 
				}
				
				# Remove existing output folder and create output directory------------
				if (file.exists(out_dir_survival)) {
					  unlink(out_dir_survival,recursive = TRUE)
				}
				dir.create(out_dir_survival, showWarnings = FALSE)
				
				## Start printing files
				progress_idx <- 1
				# loop through cancer type
				for(idx in 1:length(input_cancer_type)){
					current_cancer_type <- input_cancer_type[idx]
					# print progress
					shiny::incProgress((idx-1)/length(input_cancer_type))
					# Make outputs for each cancer type
					save_survival(
							out_dir_name  = out_dir_survival, 
							cancer_full_name = names(input_cancer_type)[input_cancer_type == current_cancer_type], 
							formatted_gene_list = plot_genes,
							labelled_data = rv[["survival"]][[current_cancer_type]],
							threshTop = survival_percent_upper(),
							threshBottom = survival_percent_lower(), 
							savePlots = FALSE,
							event = display_survival_event())
				} 
			})
			# Zip file
			shiny::withProgress(
				message = "Zipping files. This could take some time.",
				value = 0,
				{     
					shiny::incProgress(1/2)
					all_files <- list.files(out_dir_survival, full.names = TRUE)
					zip::zip(zipfile = file, files = all_files)
				}
			)
					#}, # Content function
				
				#) # Download handler
				
		},
		contentType = "application/zip" # Content 
			#) # withProgress
		) # downloadHandler

	} # function
  ) # moduleServer
} #function
