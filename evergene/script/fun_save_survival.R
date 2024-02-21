
# function that takes in dataframes required for PCA, and make all plots for all genes
	
save_survival <- function(out_dir_name, 
					cancer_full_name, 
					formatted_gene_list,
					labelled_data,
					threshTop,
					threshBottom, 
					savePlots = FALSE,
					event = "PFI"){
	#### Make directories
	dir.create(file.path(paste0(out_dir_name, "/", cancer_full_name)), showWarnings = FALSE)
	dir.create(file.path(paste0(out_dir_name, "/", cancer_full_name, "/data")), showWarnings = FALSE)

	plot_data_dir <- file.path(paste0(out_dir_name, "/", cancer_full_name, "/data/plot_data/"))
	plot_dir <- file.path(paste0(out_dir_name, "/", cancer_full_name, "/plots/"))
	dir.create(plot_data_dir, showWarnings = FALSE)
	dir.create(plot_dir, showWarnings = FALSE)

	## Create lists for saving data/plots
	stats_logrank <- list()
	stats_cox <- list()
	plot_data <- list()
	plots <- list()
	
	# Loop through genes
	for(gene in formatted_gene_list){
		gene_full_name <- names(formatted_gene_list)[formatted_gene_list == gene]
		gene_full_name <- gsub(": ", "_", gene_full_name)
		current_out <- make_survival_plot(
					cancer_full_name = cancer_full_name, 
					gene = gene, 
					gene_full_name = gene_full_name,
					labelled_data = labelled_data, 
					threshBottom = threshBottom,
					threshTop = threshTop, 
					event = event)
		stats_logrank[[gene_full_name]] = current_out[["Log_rank_stats"]]
		stats_cox[[gene_full_name]] = current_out[["Cox_stats"]]
		
		# write table for current gene
		current_out_fp <- paste0(plot_data_dir, "/", gene_full_name , ".csv")
		write.csv(current_out[["plot_data"]], current_out_fp, row.names=FALSE)
		
		# save figure for current gene
		current_out_fp <- paste0(plot_dir, "/", gene_full_name, ".pdf")
		pdf(current_out_fp ,width=6,height=6)
		#replayPlot(current_out[["surv_plot"]])
		print(current_out[["surv_plot"]])
		dev.off()

	}
	
	#### Save data tables
	out_fp_stats_logrank <-  paste0(out_dir_name, "/", cancer_full_name, "/data/stats_Log_rank.csv")
	out_fp_stats_Cox <-  paste0(out_dir_name, "/", cancer_full_name, "/data/stats_Cox_regression.csv")
	stats_logrank <- do.call(rbind, stats_logrank)
#	stats_logrank$p.adjust <- p.adjust(stats_logrank$p.value, "BH", n=nrow(stats_logrank))
	stats_cox <- do.call(rbind, stats_cox)
#	stats_cox$p.adjust <- p.adjust(stats_cox$p.value, "BH", n=nrow(stats_cox))
	
	write.csv(stats_logrank, out_fp_stats_logrank, row.names=FALSE)
	write.csv(stats_cox, out_fp_stats_Cox, row.names=FALSE)

}	
	