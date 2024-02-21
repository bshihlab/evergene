
# function that takes in dataframes required for PCA, and make all plots for all genes
	
save_pca <- function(out_dir_name, 
					cancer_type, 
					cancer_full_name, 
					scree_df, 
					pca_df, 
					top_contributions, 
					formatted_gene_list,
					add_sampleAnnot_pca,
					annotation_numeric,					
					discrete_cols,
					col_label,
					display_pcx,
					display_pcy,
					display_pca_annotation,
					savePlots = FALSE){

	shiny::withProgress(
	message = paste0("Making plots for ", cancer_full_name),
	value = 0,
	{

	#### Make directories
	dir.create(file.path(paste0(out_dir_name, "/", cancer_full_name)), showWarnings = FALSE)
	dir.create(file.path(paste0(out_dir_name, "/", cancer_full_name, "/data")), showWarnings = FALSE)
	#### Save data tables
	## Save gene contribution table
	out_fp_scree <-  paste0(out_dir_name, "/", cancer_full_name, "/data/variance_contribution.csv")
	# Because not all genes are in scree data (since low expressing genes are removed from PCA analysis)
	# Select for only overall variance explained by each PC and the %contribution from genes that have column in the scree data
	# Change the column names for the genes to gene_ensemblID so it's more readable

	formatted_gene_list_in_screedf <- formatted_gene_list[formatted_gene_list %in% colnames(scree_df)]
	if(length(formatted_gene_list_in_screedf) > 0){
		out_scree_df <- scree_df[,c("PC", formatted_gene_list_in_screedf)]
		colnames(out_scree_df)[2:ncol(out_scree_df)] <- names(formatted_gene_list_in_screedf)
	} else {
		out_scree_df <- data.frame(PC = scree_df$PC)
	}
	write.csv(out_scree_df, out_fp_scree, row.names=FALSE)

	## Save pca position/sample annotation table
	# Change the column names for the genes to gene_ensemblID so it's more readable
	out_fp_pca <-  paste0(out_dir_name, "/", cancer_full_name, "/data/pca_data.csv")
	annotate_column <- colnames(pca_df)[!colnames(pca_df) %in% formatted_gene_list]
	if(length(annotate_column)>0){
		out_pca_df <- cbind(pca_df[,annotate_column], pca_df[,formatted_gene_list])
	} else {
		out_pca_df <- pca_df[,formatted_gene_list]
	}
	colnames(out_pca_df)[(ncol(out_pca_df) - length(formatted_gene_list) + 1):ncol(out_pca_df)] <- names(formatted_gene_list)
	write.csv(out_pca_df, out_fp_pca, row.names=FALSE)	

	## Save the top contributing genes to total variation for each cancer type
	out_fp <-  paste0(out_dir_name, "/", cancer_full_name, "/data/top_gene_contributions.csv")
	current_scree_data <- scree_df[,c("PC", "contribution")]
	colnames(current_scree_data) <- c("PC", "variance.explained.by.PC") 
	contribution_df <- merge(top_contributions, current_scree_data, by="PC")
	write.csv(contribution_df, out_fp, row.names=FALSE)	
	

	#### Save plots
	if(savePlots){
		dir.create(file.path(paste0(out_dir_name, "/", cancer_full_name, "/pca_gene")), showWarnings = FALSE)
		dir.create(file.path(paste0(out_dir_name, "/", cancer_full_name, "/pc_contributions")), showWarnings = FALSE)
		dir.create(file.path(paste0(out_dir_name, "/", cancer_full_name, "/current_selection")), showWarnings = FALSE)
#		dir.create(file.path(paste0(out_dir_name, "/", cancer_full_name, "/pca_sampleAnnotation")), showWarnings = FALSE)
#		dir.create(file.path(paste0(out_dir_name, "/", cancer_full_name, "/gene_sampleAnnotation")), showWarnings = FALSE)
		
		pca_combination <- list(PC1_2 = c(1,2), PC3_4 = c(3,4), PC5_6 = c(5,6), PC7_8 = c(7,8), PC9_10 = c(9,10))
		
		progress_idx2 <- 0
		# Loop through each gene
		for(current_gene in formatted_gene_list){
			# print progress
			shiny::incProgress(progress_idx2/(length(formatted_gene_list)))
			progress_idx2 <- progress_idx2+1
			
			# organise labels
			gene_full_name <- names(formatted_gene_list)[formatted_gene_list == current_gene]
			gene_full_name <- gsub(": ", "_", gene_full_name)
			# Print scree plots
			out_plot_fp <- paste0(out_dir_name, "/",  cancer_full_name, "/pc_contributions/", gene_full_name, ".png")
			p <- make_scree_plot(current_gene, scree_df, gene_full_name, cancer_full_name)
			ggsave(out_plot_fp, width  = 4, height = 3)
			# Loop through PCs
			for(combination_name in names(pca_combination)){
				out_plot_fp <- paste0(out_dir_name, "/",  cancer_full_name, "/pca_gene/", gene_full_name , "_", combination_name, ".png")
				pcx <- pca_combination[[combination_name]][1]
				pcy <- pca_combination[[combination_name]][2]
				p <- make_pca_plot(current_gene, pcx, pcy, scree_df, pca_df, gene_full_name, cancer_full_name, col_label = col_label, discrete_cols = NA, remove_na = FALSE)
				ggsave(out_plot_fp, width = 4, height = 4)
			}
			## Plot currently selected PCs with each gene
			out_plot_fp <- paste0(out_dir_name, "/",  cancer_full_name, "/current_selection/", gene_full_name, "_PC", display_pcx, "_PC", display_pcy, ".pdf")
			p <- make_pca_plot(current_gene, display_pcx, display_pcy, scree_df, pca_df, current_gene, cancer_full_name, col_label = col_label, discrete_cols = discrete_cols, remove_na = FALSE)
			ggsave(out_plot_fp, p, width = 4, height = 4)

			## plot selected annotation
			current_annotation_label <- names(add_sampleAnnot_pca)[add_sampleAnnot_pca == display_pca_annotation]
			current_annotation_label <- gsub(": ", "_", current_annotation_label)
			# Plot each annotation against the currently selected PCs
			out_plot_fp <- paste0(out_dir_name, "/",  cancer_full_name, "/current_selection/", display_pca_annotation, "_PC", display_pcx, "_PC", display_pcy, ".pdf")
			p <- make_pca_plot(display_pca_annotation, display_pcx, display_pcy, scree_df, pca_df, current_annotation_label, cancer_full_name, col_label = "", discrete_cols = discrete_cols, remove_na = FALSE)
			ggsave(out_plot_fp, p, width = 4, height = 4)

			#Plot correlation between gene and selected annotation
			out_plot_fp <- paste0(out_dir_name, "/",  cancer_full_name, "/current_selection/", gene_full_name, "_", current_annotation_label, "_pc", display_pcy, ".pdf")
			p <- make_gene_annotation_plot (	
									pca_df = pca_df,
									scree_df = scree_df,
									input_display_gene = current_gene, 
									input_display_cancer_type = cancer_type, 
									pca_annotation = display_pca_annotation, 
									annotation_full_name = current_annotation_label, 
									gene_full_name = gene_full_name, 
									input_display_pcy = display_pcy,
									cancer_full_name = cancer_full_name, 
									remove_na = FALSE)
			ggsave(out_plot_fp, p, width = 4, height = 4)
				
		
		}
	}
	})
}	
	