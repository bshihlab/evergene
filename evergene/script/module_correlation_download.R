#### PCA: Download PDF
correlationDownloadPdf <- function(id, display_correlation_annotation1, correlation_log_gene_exp, rv) {
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
		filename = function() {"evergene_correlation.pdf"}, 
		content = function(file) {
		on.exit(removeModal())
		shiny::withProgress(
			message = "Saving - this could take some time.",
			value = 0,
			{
				# Work out the genes to plot - only plot the first pca_download_combination_restriction gene-cancer combinations
				input_cancer_type <- rv[["input_cancer_type"]]
				
				# open pdf file
				pdf(file, onefile = TRUE, width = 15, height = 10)
				progress_idx <- 1
				for(idx in 1:length(input_cancer_type)){
					current_cancer_type <- input_cancer_type[idx]
					current_cancer_full_name <- names(input_cancer_type)[idx]
					all_y_axis <- rv[["dropdown_correlation"]][[current_cancer_type]]
					all_y_axis <- all_y_axis[!all_y_axis %in% display_correlation_annotation1()]
					# maximum of 30 is selected
					if(length(all_y_axis)>30){
						all_y_axis <- all_y_axis[1:30]
					}
					# make graph 
					corr_analysed  <- make_corr_plot(
								plot_data = rv[["pca"]][[current_cancer_type]],
								x_axis = display_correlation_annotation1(), 
								y_axis = all_y_axis,
								full_names = rv[["dropdown_correlation"]][[current_cancer_type]], 
								log_gene_exp = correlation_log_gene_exp(),  
								cancer_full_name = names(rv[["input_cancer_type"]])[rv[["input_cancer_type"]] == current_cancer_type])
					print(corr_analysed[["corr_plot"]])
					# print progress
					shiny::incProgress(progress_idx/(length(input_cancer_type)))
					progress_idx <- progress_idx + 1
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
correlationDownloadZip <- function(id, display_correlation_annotation1, correlation_log_gene_exp, rv) {
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
		filename = function(){"evergene_correlation.zip"},
		content = function(file) {
			on.exit(removeModal())
			shiny::withProgress(
			message = "Saving - this could take some time.",
			value = 0,
			{
				# Work out the genes to plot - only plot the first pca_download_combination_restriction gene-cancer combinations
				input_cancer_type <- rv[["input_cancer_type"]]
				
				
				# Remove existing output folder and create output directory------------
				if (file.exists(out_dir_correlation)) {
					  unlink(out_dir_correlation,recursive = TRUE)
				}
				dir.create(out_dir_correlation, showWarnings = FALSE)
				
				## Start printing files
				progress_idx <- 1
				# loop through cancer type
				for(idx in 1:length(input_cancer_type)){
					current_cancer_type <- input_cancer_type[idx]
					current_cancer_name <- names(input_cancer_type)[idx]
					all_y_axis <- rv[["dropdown_correlation"]][[current_cancer_type]]
					# maximum of 30 is selected
					if(length(all_y_axis)>30){
						all_y_axis <- all_y_axis[1:30]
					}
					current_outdir <- paste0(out_dir_correlation, "/", current_cancer_name)
					dir.create(current_outdir, showWarnings = FALSE)
					
					# Make outputs for each cancer type
					corr_analysed  <- make_corr_plot(
								plot_data = rv[["pca"]][[current_cancer_type]],
								x_axis = display_correlation_annotation1(), 
								y_axis = all_y_axis,
								full_names = rv[["dropdown_correlation"]][[current_cancer_type]], 
								log_gene_exp = correlation_log_gene_exp(),  
								cancer_full_name = names(rv[["input_cancer_type"]])[rv[["input_cancer_type"]] == current_cancer_type])
					ggsave(paste0(current_outdir, "/plot.pdf"), plot=corr_analysed[["corr_plot"]], width = 15, height = 10)
					write.csv(corr_analysed[["correlation"]], paste0(current_outdir, "/correlation.csv"), row.names=FALSE)
					write.csv(corr_analysed[["plot_df"]], paste0(current_outdir, "/plot_data.csv"), row.names=FALSE)
					
					# print progress
					shiny::incProgress(progress_idx/(length(input_cancer_type)))
					progress_idx <- progress_idx + 1
				} 	
			})
			# Zip file
			shiny::withProgress(
				message = "Zipping files. This could take some time.",
				value = 0,
				{     
					shiny::incProgress(1/2)
					all_files <- list.files(out_dir_correlation, full.names = TRUE)
					zip::zip(zipfile = file, files = all_files)
				}
			)			
		},
		contentType = "application/zip" # Content 
		) # downloadHandler
	} # function
  ) # moduleServer
} #function

