# takes in data for plotting PCA
make_corr_plot <- function( plot_data,
							x_axis,
							y_axis,
							full_names,
							log_gene_exp,
							cancer_full_name,
							show_corrTest = TRUE,
							asterisk_pval = 10^(-20)){
	# Organise dataframe for plotting

	## organise if x-axis needs to be converted
	plot_data <- as.data.frame(plot_data)
	colnames(plot_data)[1] <- "sample_id"
	xlabel <- names(full_names)[as.character(full_names) == x_axis]
	corr_x <- as.numeric(plot_data[,x_axis])

	# log values for correlation
	if(log_gene_exp & length(grep(": ENSG0", xlabel)) > 0){
		# add 0.1 TPM of the smallest non-zero value before log (to get around log 0)
		corr_x <- corr_x + 0.1
		corr_x <- log10(corr_x)
		xlabel <- paste0( xlabel, "(log10TPM)")
	} else if (length(grep(": ENSG0", xlabel)) > 0){
		xlabel <- paste0( xlabel, "(TPM)")
	}
	
	## Organise each y-axis annotation
	# there may be more than 1 y-axis selected
	plot_df_list <- list()
	y_lab_factor <- vector() # for ordering y_lab
	for (current_y in y_axis){
		current_y_label <- names(full_names)[as.character(full_names) == current_y]
		# remove the ensembl gene ID for plots
		current_df <- data.frame(sample_id = plot_data$sample_id, x=corr_x, y=plot_data[,current_y], y_id = current_y, y_lab_full = current_y_label)
		current_df$y <- as.numeric(current_df$y)
		# log10 y if it is a gene and log_gene_exp is on
		if(log_gene_exp & length(grep(": ENSG0", current_y_label)) > 0){
			# add 0.1 TPM of the smallest non-zero value before log (to get around log 0)
			current_df$y <- current_df$y + 0.1
			current_df$y <- log10(current_df$y)
			current_y_label <- strsplit(current_y_label, ": ENS")[[1]][1]
			current_y_label <- paste0( current_y_label, "(log10TPM)")
		} else if (length(grep(": ENSG0", current_y_label)) > 0){
			current_y_label <- strsplit(current_y_label, ": ENS")[[1]][1]
			current_y_label <- paste0( current_y_label, "(TPM)")
		}
		# store the editted y-label
		current_df$y_lab <- current_y_label
		# storing input y-axis order 
		y_lab_factor <- c(y_lab_factor, current_y_label)
	
		plot_df_list[[current_y]] = current_df	
	
	}
	plot_df <- do.call(rbind, plot_df_list)
	plot_df$y_lab <- factor(plot_df$y_lab, levels=y_lab_factor)
	
	## remove any datapoints with missing x/y
	plot_df <- plot_df[!is.na(plot_df$x) & !is.na(plot_df$y),]

	## Labels for the title
	plot_title <- paste0(cancer_full_name, "\n")
	
	## Calculate correlation
	if(show_corrTest){
		# get a vector for the correlation test results
		corr_output <- lapply(names(plot_df_list), function(x){
			current_df <- plot_df_list[[x]]
			pearson_r <- cor(current_df[,"x"], current_df[,"y"], use="pairwise.complete.obs")
			pearson_r_pretty <- prettyNum(pearson_r, digits = 2)
			pearson_p <- cor.test(current_df[,"x"], current_df[,"y"], na.action = "na.omit")
			pearson_p_pretty <- prettyNum(pearson_p$p.value, digits = 2)
			plot_text <- ifelse(pearson_p$p.value < asterisk_pval, paste0("(r=", pearson_r_pretty, ", p.val=", pearson_p_pretty, ")*"), paste0("(r=", pearson_r_pretty, ", p.val=", pearson_p_pretty, ")"))
			return(data.frame(r = pearson_r, p_val = pearson_p$p.value, corr_results = plot_text))
		})
		corr_df <- do.call(rbind, corr_output)
		corr_df$y_id <- names(plot_df_list)
		plot_df <- merge(plot_df, corr_df, by="y_id")


		# add correlation results
		plot_df$y_lab_new <- paste0(plot_df$y_lab, "\n", plot_df$corr_results)
		plot_df <- plot_df[order(plot_df$y_lab),]
		plot_df$y_lab <- factor(plot_df$y_lab_new, levels=unique(plot_df$y_lab_new))
		# clean up output table
		plot_df$x_lab <- xlabel
		plot_df$y_lab_new <- NULL
	}


	## Plot
	out_plot <- ggplot(data = plot_df, aes(x = x, y = y, label=sample_id)) + 
		geom_point(size=0.8) +
	  theme_bw() +  
	  facet_wrap(~y_lab, scales = "free_y") +
	  labs(x = xlabel, 
		  y = "[y-axis and correlation test indicated in graph titles]",
		  title=plot_title) +
	  theme(plot.margin=unit(c(8,2,0,0),"mm"), 
		text=element_text(size=10),
		plot.title = element_text(size=12, margin=margin(0,0,0,0)),
		panel.grid=element_blank(),
		strip.text.x = element_text(margin = margin(0.2,0,0.2,0, "cm")))

	## tidy up output dataframe
	correlation_df_out <- plot_df
	correlation_df_out$x=NULL
	correlation_df_out$y=NULL
	correlation_df_out <- aggregate(data = correlation_df_out, y_id ~ x_lab + y_lab + y_lab_full + r + p_val, length)
	colnames(correlation_df_out)[ncol(correlation_df_out)] <- "sample_count"
	
	# reorder the rows
	correlation_df_out <- correlation_df_out[order(correlation_df_out$y_lab),]
	
	return(list(corr_plot = out_plot, correlation = correlation_df_out, plot_df = plot_data))
}
