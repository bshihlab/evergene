# takes in data for plotting PCA
make_gene_annotation_plot <- function(
									pca_df,
									scree_df, 
									input_display_gene,
									input_display_cancer_type,
									pca_annotation,
									annotation_full_name,
									gene_full_name, 
									cancer_full_name, 
									input_display_pcy ,
									remove_na = FALSE){
	# Get data
	# use gene name to label scree plot fill legend if there is a gene name (if not, use gene_id)
	# make correlation plot
	plot_df <- pca_df
	colnames(plot_df)[1] <- "sample_id"
	plot_df[[input_display_gene]] <- plot_df[[input_display_gene]] + 0.1 # add 0.1 to allow log scaling
	# Check if current annotation should be numeric
	if(pca_annotation %in% annotation_numeric){
		plot_df[[pca_annotation]] <- as.numeric(plot_df[[pca_annotation]])
	}	
	is_numeric <- is.numeric(plot_df[[pca_annotation]])

	
	# if input_display_pcy is NA, use the PC with the highest contribution for the input gene
	if(is.na(input_display_pcy)){
		if(input_display_gene %in% colnames(scree_df) ){
			input_display_pcy <- scree_df[[input_display_gene]]
			input_display_pcy <- scree_df$PC[which.max(input_display_pcy)]
		} else {
			input_display_pcy <- "PC1"
		}
	} else {
		input_display_pcy <- paste0("PC",input_display_pcy)
	}

	# Labels for the plot title
	plot_title <- paste0(gene_full_name, "\n", annotation_full_name, "\n", cancer_full_name, " PC [colour]\n")

	# Check if jitter should be used for current annotation
	jitter_condition <- ifelse(pca_annotation %in% annotation_jitter, TRUE, FALSE)
	
	# use pre-processed data to make dotplot
	corr_plot <- ggplot(plot_df, aes(x =.data[[input_display_gene]], y=.data[[pca_annotation]], colour=.data[[input_display_pcy]], label=sample_id)) + 
				theme_bw() + xlab(paste0(gene_full_name)) + 
				theme(plot.margin=unit(c(14.5,0,0,0),"mm"), 
					text=element_text(size=10),
					plot.title = element_text(size=12, margin=margin(0,0,0,0)),
					panel.grid=element_blank()) + scale_x_log10()
					
					
	# Add points depending on if jittering is required
	if(jitter_condition){
		corr_plot <- corr_plot+ geom_point(position = position_jitter(width = 0.1, height = 0.2)) 
	} else {
		corr_plot <- corr_plot+ geom_point()
	}
	# Add point colour
	if(is_numeric){
		corr_plot <- corr_plot + scale_colour_viridis(name = input_display_pcy, option = "A", direction = -1, end=0.93) + 
						ylab(annotation_full_name) + 
						theme(plot.margin=unit(c(14.5,0,0,0),"mm"), 
							text=element_text(size=10),
							plot.title = element_text(size=12, margin=margin(0,0,0,0))) +
						labs(title = plot_title)
		ggplotly(corr_plot)
	} else {
		plot_cols <- discrete_cols[1:length(unique(plot_df[[pca_annotation]]))]
		corr_plot <- corr_plot + scale_colour_manual(name = "annotation", values = plot_cols) + 
					ylab("") + 
					stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", colour="#000000", width=0.4, linewidth=0.3)	+ 			
					stat_summary(fun = function(x)quantile(x)[2], fun.min = function(x)quantile(x)[2], fun.max = function(x)quantile(x)[2], geom = "crossbar", colour="#000000", width=0.4, linewidth=0.1)	+		
					stat_summary(fun = function(x)quantile(x)[4], fun.min = function(x)quantile(x)[4], fun.max = function(x)quantile(x)[4], geom = "crossbar", colour="#000000", width=0.4, linewidth=0.1)	+
					theme(plot.margin=unit(c(14.5,0,0,0),"mm"), 
						text=element_text(size=10),
						plot.title = element_text(size=12, margin=margin(0,0,0,0)))	+
					scale_colour_viridis(name = input_display_pcy, option = "A", direction = -1, end=0.93) +
					labs(title = plot_title)
	}	

	return(corr_plot)
}
