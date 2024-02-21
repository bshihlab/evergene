#### Make barplots showing number of selected projects, detected genes, and invalid genes

make_gene_project_count_plot <- function(gene_string, projects, combination_cutoff = pca_plot_combination_restriction){
	# Organise plotting dataframe
	formatting_output <- input_format(gene_string, projects)
	genes_found_num <- length(formatting_output[["formatted_gene_list_named"]])
	genes_notFound_num <- length(formatting_output[["rejected_gene_list"]])
	project_num <- length(projects)
	
	
	## Plot:
	# - selected number of cancer projects
	# - selected number of genes (recognised)
	# - selected number of genes (not recognised)
	
	plot_data <- data.frame(input_data=c("Cancer type", "Genes found", "Genes not found"), 
							count=c(project_num, genes_found_num, genes_notFound_num),
							label = c("Input cancer types", 
										"Input genes recognised",
										"Input genes not recognised"))

	# Check the total number of gene-project combinations, and whether there are any invalid genes
	# Colour code the messages
	gene_project_combinations <- genes_found_num*project_num
	if(gene_project_combinations <= combination_cutoff & genes_notFound_num == 0){
		output_message <- paste0(gene_project_combinations, " gene-cancer combinations found\n Ready for ANALYSE")
		output_message_col <- "#5FA052" # muted green
	} else if(gene_project_combinations <= combination_cutoff & genes_notFound_num > 0 & genes_found_num > 0 ) {
		output_message <- paste0(gene_project_combinations, " gene-cancer combinations found\n Ready for ANALYSE\n Note: ", genes_notFound_num, " invalid genes")
		output_message_col <- "#EB9C5C" # muted orange	
	} else {
		output_message <- paste0(gene_project_combinations, " gene-cancer combinations found\n Please keep it between 1-", combination_cutoff)
		output_message_col <- "#CB4C4E"	# muted red
	}
	
	# make plot

	p <- ggplot(plot_data,aes(x = input_data, y = count, label=label)) +
		geom_bar(stat="identity") + theme_bw() +
		ggtitle(output_message) + xlab("") + ylab("Count") + 
		theme(
		  legend.position="none",
		  panel.background=element_blank(),
		  panel.grid.minor=element_blank(),
		  plot.background=element_blank()) +
		  scale_y_continuous(breaks=pretty_breaks()) +
#		  theme(panel.border = element_blank()) + #, axis.line = element_line()) +
		  theme(plot.margin=unit(c(8,0,0,0),"mm"), 
			text=element_text(size=10), 
			plot.title = element_text(size=10, margin=margin(0,0,0,0), colour=output_message_col)) + 
		  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	return(p)
}