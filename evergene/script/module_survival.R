## Survival main plot
survivalPlot1 <- function(id, display_cancer_type, display_gene, survival_percent_lower, survival_percent_upper, display_survival_event, rv) {
  moduleServer( id,
	function(input, output, session) {
		output$showPlot <- renderPlot({
			validate(need(display_cancer_type(), message = "cancer type missing"), 
				   need(display_gene(), message = "gene missing"),
				   need(survival_percent_lower(), message = "lower missing"),
				   need(survival_percent_upper(), message = "upper missing"),
				   need(display_survival_event(), message = ""),
				   need(survival_percent_lower() <= survival_percent_upper(), message = "Lower threshold must be less than upper threshold")) # Validate needs
				   
			# use pre-processed data to make PCA plot
			out_plot <- make_survival_plot(
						cancer_full_name = names(rv[["input_cancer_type"]])[rv[["input_cancer_type"]] == display_cancer_type()],
						gene = display_gene(),
						gene_full_name = names(rv[["formatted_gene_list"]])[rv[["formatted_gene_list"]] == display_gene()],
						labelled_data = rv[["survival"]][[display_cancer_type()]],
						plot_cols = discrete_cols_surv,
						threshBottom = survival_percent_lower(),
						threshTop = survival_percent_upper(),
						event = display_survival_event())
			out_plot[["surv_plot"]]
		})
	}
  )
}


survivalPlot2 <- function(id, display_cancer_type, display_gene, survival_percent_lower, survival_percent_upper, 
						display_survival_event, display_survival_annotation, survival_annotation_na_rm, 
						display_pcy, survival_annotation_figLegend, rv) {
  moduleServer( id,
	function(input, output, session) {
		output$showPlot <- renderPlotly({
			validate(need(display_cancer_type(), message = "cancer type missing"), 
				   need(display_gene(), message = "gene missing"),
				   need(survival_percent_lower(), message = "lower missing"),
				   need(survival_percent_upper(), message = "upper missing"),
				   need(display_survival_annotation(), message=FALSE),
				   need(display_survival_event(), message = ""),
				   need(survival_percent_lower() <= survival_percent_upper(), message = "Lower threshold must be less than upper threshold")) # Validate needs
				   
			current_pcy <- paste0("PC", display_pcy())
			current_gene_fullname <- names(rv[["formatted_gene_list"]])[rv[["formatted_gene_list"]] == display_gene()]
			plot_df <- rv[["survival"]][[display_cancer_type()]]
			# get exprs grouping (high/low expressing groups)
			plot_df <- survival_grouping(
					barcode = plot_df$barcode,
					exprs = plot_df[[display_gene()]], 
					survival_status  = as.numeric(plot_df[,display_survival_event()]), 
					survival_time  = as.numeric(plot_df[, paste0(display_survival_event(), ".time")]),
					threshBottom = survival_percent_lower() , 
					threshTop = survival_percent_upper())
			plot_df <- merge(plot_df[,c("barcode", "exprs_grp")], rv[["pca"]][[display_cancer_type()]], by="barcode")
			plot_df$sample_id <- plot_df$barcode
			print(head(plot_df))
			# Check if the values are numeric
			if(display_survival_annotation() %in% annotation_numeric){
				plot_df[[display_survival_annotation()]] <- as.numeric(plot_df[[display_survival_annotation()]])
			}

			# remove na if check box ticked
			if(survival_annotation_na_rm()){
				plot_df <- plot_df[!is.na(plot_df[[display_survival_annotation()]]),]
			}
			
			is_numeric <- is.numeric(plot_df[[display_survival_annotation()]][!is.na(plot_df[[display_survival_annotation()]])])
			
			# make plot1 PC-y vs gene exprs coloured by sample annotation 
			p1 <- ggplot(data = plot_df, 
				aes(x= .data[[display_gene()]], y =.data[[current_pcy]], 
					colour= .data[[display_survival_annotation()]]), label=sample_id) +
				geom_point(size=1) + theme_bw() +
				xlab(current_gene_fullname) +
				theme( panel.grid=element_blank())
			# change colour scales according to if it is numeric
			if(is_numeric){
				p1 <- p1 + scale_colour_viridis(option = "A", direction = -1, end=0.93)
			} else {
				p1 <- p1 + scale_colour_manual(values = discrete_cols[1:length(unique(plot_df[[display_survival_annotation()]]))]) 
			}
			

			# make plot2 PC-y vs gene exprs coloured by gene expression groups 
			p2 <- ggplot(data = plot_df, 
				aes(x= .data[[display_gene()]], y =.data[[current_pcy]], 
					colour= .data[["exprs_grp"]]), label=sample_id) +
				geom_point(size=0.8) + theme_bw() + 
				xlab(current_gene_fullname) +
				guides(colour = guide_legend(override.aes = list(size=5), title=""))  +
				theme( panel.grid=element_blank()) +
				scale_colour_manual(values = discrete_cols_exprs[unique(plot_df$exprs_grp)])
			if(min(plot_df[[display_gene()]])>0){	
				p1 <- p1 + scale_x_log10() 
				p2 <- p2 + scale_x_log10() 
			}
			if(survival_annotation_figLegend()){
				subplot(ggplotly(p1) %>% layout(showlegend = F), ggplotly(p2) %>% layout(showlegend = F), nrows=2, shareX = TRUE, shareY=TRUE )
			} else {
				subplot(ggplotly(p1), ggplotly(p2), nrows=2, shareX = TRUE, shareY=TRUE )
			}
		})
	}
  )
}
