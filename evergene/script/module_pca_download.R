#### PCA: Download PDF
pcaDownloadPdf <- function(id, rv) {
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
		filename = function() {"evergene_pca.pdf"}, 
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
				pdf(file, onefile = TRUE, width = 4.5, height = 4.5)
				progress_idx <- 1
				for(idx in 1:length(input_cancer_type)){
					current_cancer_type <- input_cancer_type[idx]
					current_cancer_full_name <- names(input_cancer_type)[idx]
					# print progress
					current_plots <- make_pca_plot_multiGene(cancer_full_name = current_cancer_full_name,
															scree_df = rv[["scree"]][[current_cancer_type]], 
															pca_df = rv[["pca"]][[current_cancer_type]], 
															formatted_gene_list = plot_genes)
					for(plot_idx in 1:length(current_plots)){
						print(current_plots[[plot_idx]])
						  shiny::incProgress(progress_idx/(length(input_cancer_type)*length(plot_genes)))
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
pcaDownloadZip <- function(id, display_pcx, display_pcy, display_pca_annotation, rv) {
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
		filename = function(){"evergene_pca.zip"},
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
				  plot_genes <- rv[["formatted_gene_list"]] 
				}
				
				# Remove existing output folder and create output directory------------
				if (file.exists(out_dir_pca)) {
					  unlink(out_dir_pca,recursive = TRUE)
				}
				dir.create(out_dir_pca, showWarnings = FALSE)
				
				## Start printing files
				progress_idx <- 1
				# loop through cancer type
				for(idx in 1:length(input_cancer_type)){
					current_cancer_type <- input_cancer_type[idx]
					# print progress
					shiny::incProgress((idx-1)/length(input_cancer_type))
					# Make outputs for each cancer type
					save_pca(
						out_dir_name = out_dir_pca, 
						cancer_type = current_cancer_type, 
						cancer_full_name = names(input_cancer_type)[input_cancer_type == current_cancer_type], 
						scree_df = rv[["scree"]][[current_cancer_type]], 
						pca_df = rv[["pca"]][[current_cancer_type]],
						top_contributions = rv[["top_contributions"]][[current_cancer_type]], 
						formatted_gene_list = plot_genes, 
						add_sampleAnnot_pca = add_sampleAnnot_pca,
						annotation_numeric = annotation_numeric,					
						discrete_cols = discrete_cols,
						col_label = rv[["exprs_unit"]],
						display_pcx = display_pcx(),
						display_pcy = display_pcy(),
						display_pca_annotation = display_pca_annotation(),
						savePlots = TRUE)
				} 
			})
			# Zip file
			shiny::withProgress(
				message = "Zipping files. This could take some time.",
				value = 0,
				{     
					shiny::incProgress(1/2)
					all_files <- list.files(out_dir_pca, full.names = TRUE)
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

