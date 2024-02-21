#### use pre-processed data to make PCA/scree plot on all combinations
# make pc pairs so it's PC1-PC2, PC3-PC4... so on
make_pca_plot_multiGene <- function(cancer_full_name, scree_df, pca_df, formatted_gene_list){
	pca_combination <- list(PC1_2 = c(1,2), PC3_4 = c(3,4), PC5_6 = c(5,6), PC7_8 = c(7,8), PC9_10 = c(9,10))
	out_plot_list <- list()
	plot_idx <- 1
	for(current_gene in formatted_gene_list){
		# organise labels
		gene_full_name <- names(formatted_gene_list)[formatted_gene_list == current_gene]
		gene_full_name <- gsub(": ", "_", gene_full_name)
		# Print scree plots
		out_plot_list[[plot_idx]] <- make_scree_plot(current_gene, scree_df, gene_full_name, cancer_full_name)
		plot_idx <- plot_idx + 1

		# Loop through PCs
		for(combination_name in names(pca_combination)){
			pcx <- pca_combination[[combination_name]][1]
			pcy <- pca_combination[[combination_name]][2]
			p <- make_pca_plot(current_gene, pcx, pcy, scree_df, pca_df, gene_full_name, cancer_full_name, col_label = "TPM", discrete_cols = NA, remove_na = FALSE)
			out_plot_list[[plot_idx]] <- p
			plot_idx <- plot_idx + 1
		}
	}
	return(out_plot_list)
}

		