# use pre-processed data to make PCA plot
make_geneContribution_plot <- function(gene, scree_df){
	# Plotting for expressed genes
	if(gene %in% colnames(scree_df)){
		plot_df <- scree_df[,c("PC", gene)]
		plot_df$direction <- sign(plot_df[,gene])
		#plot_df <- plot_df[as.numeric(plot_df$PC) <=num_pc,]
		plot_df$direction <- ifelse(plot_df$direction == 1, "positive", "negative")
		plot_df$gene <- abs(plot_df[,gene])
		tile_plot_gene <- ggplot(data = plot_df, mapping = aes(x = PC, y = 1, label = PC, fill = gene)) +	
		  geom_tile() +
		  geom_text(aes(colour = direction), size=3) +
		  scale_x_discrete(limit = paste0("PC", 1:length(unique(plot_df$PC)))) +
		  theme_minimal() + 
		  scale_fill_viridis(option = "magma", direction = -1) + xlab("") + ylab("") + 
		  scale_colour_manual(values = c("#FFFFFF", "#FFFFFF"))
	  
  	# Plotting for genes with almost no expression
	} else {
		current_scree_data <- scree_df[,c("PC", "contribution")]
		colnames(current_scree_data) <- c("PC", "contribution") 
	  
		tile_plot_gene <- ggplot(data = current_scree_data, mapping = aes(x = PC, y = 1, label = PC, )) +
		  geom_tile(fill = "grey")  +  geom_text(colour = "#ffffff", size=3) +
		  scale_x_discrete(limit = paste0("PC", 1:length(unique(plot_df$PC)))) +
		  theme_minimal() + xlab("") + ylab("")
	} # Else   
	
	tile_plot_gene <- tile_plot_gene + 
		theme( legend.position = "none" )+ 
		theme(axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                plot.margin=unit(c(0,0,0,0),"mm"),
                plot.background = element_blank(),
                panel.grid=element_blank(),
                panel.border=element_blank(),
                axis.ticks=element_blank()) +
				labs(x=NULL, y=NULL)

	return(tile_plot_gene)
}


