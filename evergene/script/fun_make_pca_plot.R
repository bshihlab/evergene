# takes in data for plotting PCA
make_pca_plot <- function(pca_annotation, pcx, pcy, scree_df, pca_df, gene_full_name, 
					cancer_full_name, col_label , discrete_cols = NA, remove_na = FALSE){
	
	# Organise dataframe for plotting
	## plot PCA
	pcx <- paste0("PC", pcx)
	pcy <- paste0("PC", pcy)
	pcx_contribution <- round(scree_df$contribution[scree_df$PC == pcx], 1) # labeling axis
	pcy_contribution <- round(scree_df$contribution[scree_df$PC == pcy], 1)  # labeling axis
	
	# main plot data
	current_pca_plot_df <- data.frame(
		sample_id = pca_df[[colnames(pca_df)[1]]],
		val = pca_df[[pca_annotation]],
		pcx = pca_df[[pcx]],
		pcy = pca_df[[pcy]])

	# Labels for the title
	plot_title <- paste0(cancer_full_name, " PCs\n", gene_full_name, " [colour]")
	
	if(pca_annotation %in% annotation_numeric){
		current_pca_plot_df$val <- as.numeric(current_pca_plot_df$val)
	}
	
	# Check if the values are numeric
	is_numeric <- is.numeric(current_pca_plot_df$val[!is.na(current_pca_plot_df$val)])
	
	# remove rows with na
	if(remove_na){
		current_pca_plot_df <- current_pca_plot_df[!is.na(current_pca_plot_df$val),]
	}

	# change column names
	# Only log expression if it's tcga data, as it is not known what data type the user input is
	if((col_label== "TPM") & (cancer_full_name != "User input")){
		# if the fill colour is a gene, add 0.1 to the TPM so it can be log transformed
		current_pca_plot_df$val <- current_pca_plot_df$val + 0.1
		colnames(current_pca_plot_df) <- c("sample_id", paste0(pca_annotation, ".expression.tpm"), "pcx", "pcy")
	} else {
		colnames(current_pca_plot_df) <- c("sample_id", pca_annotation, "pcx", "pcy")
	}
	plot_name <- colnames(current_pca_plot_df)[2]

	# Plot
	pca_plot <- ggplot(data = current_pca_plot_df, aes(x = pcx, y = pcy, label=sample_id)) + 
	  labs(title = plot_title, colour = col_label) + 
	  theme_bw() + 
	  xlab(paste0(pcx, " (", pcx_contribution, "%)")) + 
	  ylab(paste0(pcy, " (", pcy_contribution, "%)")) +
	  theme(plot.margin=unit(c(7,0,0,0),"mm"), 
		text=element_text(size=10),
		plot.title = element_text(size=12, margin=margin(0,0,0,0)),
		panel.grid=element_blank())

	# account for non-numeric (discrete colouring)
	if(is_numeric){
		if((col_label== "TPM") & (cancer_full_name != "User input")){
			my_breaks = round(exp(seq(log(1), log(max(current_pca_plot_df[[plot_name]])), length=6)), 0)
			pca_plot <- pca_plot + geom_point(aes(colour = .data[[plot_name]])) + scale_colour_viridis(option = "A", trans = "log", breaks = my_breaks, labels = my_breaks, direction = -1, end=0.95)
		} else {
			pca_plot <- pca_plot + geom_point(aes(colour = .data[[plot_name]])) + scale_colour_viridis(option = "A", direction = -1, end=0.95)
		}
	} else {
		discrete_cols <- discrete_cols[1:length(unique(current_pca_plot_df[[plot_name]]))]
		pca_plot <- pca_plot + geom_point(aes(colour = .data[[plot_name]])) + scale_colour_manual(values = discrete_cols)
	}
  

	return(pca_plot)
}
