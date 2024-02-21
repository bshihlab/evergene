## Correlation plot
correlationPlot1 <- function(id, display_cancer_type, display_correlation_annotation1, display_correlation_annotation2, corr_plot_all_options, correlation_log_gene_exp, rv) {
  moduleServer( id,
	function(input, output, session) {
		output$showPlot <- renderPlotly({
			validate(need(display_cancer_type(), message = FALSE), 
				   need(display_correlation_annotation1(), message = "annotation1_missing"),
				   need(display_correlation_annotation2(), message = "annotation2_missing"),
				   need(!rv$computingInProgress, message = "computing in progress"))

			# use pre-processed data to make PCA plot
			corr_analysed <- make_corr_plot(
						plot_data = rv[["pca"]][[display_cancer_type()]],
						x_axis = display_correlation_annotation1(), 
						y_axis = display_correlation_annotation2(),
						full_names = rv[["dropdown_correlation"]][[display_cancer_type()]], 
						log_gene_exp = correlation_log_gene_exp(),  
						cancer_full_name = names(rv[["input_cancer_type"]])[rv[["input_cancer_type"]] == display_cancer_type()])
			ggplotly(corr_analysed[["corr_plot"]]) 
		})
	}
  )
}

