# takes in data for plotting PCA
make_topGeneContribution_plot <- function( 
						contribution_df,
						scree_df, 
						top_gene_contribution){
	# Get data
	contribution_df <- contribution_df[contribution_df$gene_rank <= top_gene_contribution,]
	current_scree_data <- scree_df[,c("PC", "contribution")]
	colnames(current_scree_data) <- c("PC", "variance.explained.by.PC") 
	contribution_df <- merge(contribution_df, current_scree_data, by="PC")
	contribution_df$PC <- factor(contribution_df$PC, levels=paste0("PC", 1:length(unique(contribution_df$PC))))
	contribution_df$contribution.whole.dataset <- contribution_df$variance.explained.by.PC * abs(contribution_df$contribution)/100
	
	# reverse factor so the top rank comes in the first row
	contribution_df$gene_rank <- factor(contribution_df$gene_rank, levels=rev(1:top_gene_contribution))
	
	# use gene name to label scree plot fill legend if there is a gene name (if not, use gene_id)
	# make tile plot
	plot_df <- contribution_df[,c("PC", "gene_rank", "gene", "contribution.whole.dataset",  "contribution")]
	plot_df$direction <- sign(plot_df$contribution)
	plot_df$direction <- ifelse(plot_df$direction == 1, "positive", "negative")
	colour_bins <- seq(0,max(plot_df$contribution.whole.dataset), length.out=10)
	out_plot <- ggplot(plot_df, aes(x=PC, y = gene_rank, fill=contribution.whole.dataset, label=gene)) +
				geom_tile() + ylab("Rank") + xlab("") +
				geom_text(size = 2.5, aes(colour = direction )) + 
				scale_fill_viridis(option = "A", direction = -1, end=0.93, guide = "none") + theme_bw() +
				theme(plot.margin=unit(c(0,0,0,0),"mm"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
				scale_colour_manual(values=c("#FFFFFF", "#ffffff"))+ 
				theme(legend.position = "bottom")

	return(out_plot)
}


# takes in data for plotting PCA
make_topGeneContribution_plot_archive <- function( 
						contribution_df,
						scree_df, 
						top_gene_contribution){
	# Get data
	contribution_df <- contribution_df[contribution_df$gene_rank <= top_gene_contribution,]
	current_scree_data <- scree_df[,c("PC", "contribution")]
	colnames(current_scree_data) <- c("PC", "variance.explained.by.PC") 
	contribution_df <- merge(contribution_df, current_scree_data, by="PC")
	contribution_df$PC <- factor(contribution_df$PC, levels=paste0("PC", 1:length(unique(contribution_df$PC))))
	contribution_df$contribution.whole.dataset <- contribution_df$variance.explained.by.PC * abs(contribution_df$contribution)
	
	# reverse factor so the top rank comes in the first row
	contribution_df$gene_rank <- factor(contribution_df$gene_rank, levels=rev(1:top_gene_contribution))
	
	# use gene name to label scree plot fill legend if there is a gene name (if not, use gene_id)
	# make tile plot
	plot_df <- contribution_df[,c("PC", "gene_rank", "gene", "contribution.whole.dataset",  "contribution")]
	plot_df$direction <- sign(plot_df$contribution)
	plot_df$direction <- ifelse(plot_df$direction == 1, "positive", "negative")
	out_plot <- ggplot(plot_df, aes(x=PC, y = gene_rank, fill=contribution.whole.dataset , label=gene)) +
				geom_tile() + ylab("Rank") + xlab("") +
				geom_text(size = 2.5, aes(colour = direction )) + 
				scale_fill_viridis(option = "A", direction = -1, end=0.93) + theme_bw() +
				theme(plot.margin=unit(c(0,0,0,0),"mm"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
				scale_x_discrete(position = "top") +
				scale_colour_manual(values=c("#FFFFFF", "#ffffff"))
	return(out_plot)
}
