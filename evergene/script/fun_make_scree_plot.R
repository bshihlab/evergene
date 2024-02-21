# use pre-processed data to make PCA plot
make_scree_plot <- function(gene, scree_df, gene_full_name, cancer_full_name){
	# Plotting for expressed genes
	if(gene %in% colnames(scree_df)){
		title_label <- paste0("Scree plot\n", cancer_full_name,"\n", gene_full_name, "\n   contribution to each PC [colour]")
		current_scree_data <- scree_df[,c("PC", "contribution", gene)]
		colnames(current_scree_data) <- c("PC", "variance.explained.by.PC", "gene.contribution") 
		current_scree_data$PC <- gsub("PC", "", current_scree_data$PC)
		current_scree_data$PC  <- factor(current_scree_data$PC , levels=1:length(unique(current_scree_data$PC)))
		current_scree_data$direction <- sign(current_scree_data$gene.contribution)
		current_scree_data$gene.contribution  <- abs(current_scree_data$gene.contribution)
	  
		scree_plot <- ggplot(data = current_scree_data, mapping = aes(x = PC, y = variance.explained.by.PC, fill = gene.contribution)) +
			geom_bar(stat = "identity", colour="black", linewidth=0.1)  + theme_bw() +
			scale_fill_viridis(name = "Contribution(%)",  option = "magma", direction = -1) +
			theme(plot.margin=unit(c(28,0,0,0),"mm"), 
				text=element_text(size=10),
				plot.title = element_text(size=12, margin=margin(0,0,0,0)),
				panel.grid=element_blank())+
				labs(title = title_label) +
				xlab("Principal Component (PC)") + ylab("Variance Explained by PC (%)") 
		# Plotting for genes with almost no expression
	} else {
		title_label <- paste0("Scree plot\n", cancer_full_name)
		current_scree_data <- scree_df[,c("PC", "contribution")]
		colnames(current_scree_data) <- c("PC", "variance.explained.by.PC") 
		
		scree_plot <- ggplot(data = current_scree_data, mapping = aes(x = PC, y = variance.explained.by.PC)) +
		  geom_bar(stat = "identity", colour="black", linewidth=0.1, fill="grey") + theme_bw() +
		  scale_x_discrete(limit = paste0("PC", length(unique(current_scree_data$PC)))) +
		  theme(plot.margin=unit(c(10,0,0,0),"mm"), 
				text=element_text(size=10),
				plot.title = element_text(size=12, margin=margin(0,0,0,0)),
				panel.grid=element_blank())+
				labs(title = title_label) +
				xlab("Principal Component (PC)") + ylab("Variance Explained by PC (%)") 
	} # Else   
	
 
	return(scree_plot)
}