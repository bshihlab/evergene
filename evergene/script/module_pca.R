### PCA gene contribution tile plot at the top
pcaGeneContributionServer <- function(id, display_cancer_type, display_gene, rv) {
	moduleServer( id,
	function(input, output, session) {
	  output$showPlot <- renderPlotly({
			# Validate needs
			validate(need(display_cancer_type(), message = FALSE), need(display_gene(), message =  FALSE))
			# make tile plot on the top, shared between PCA and survival tabs
			out_plot <- make_geneContribution_plot(
				gene = display_gene(),
				scree_df = rv[["scree"]][[display_cancer_type()]])		
			ggplotly(out_plot)				
		})
	})
}

## PCA plot1 (colour by selected gene)
pcaPlot1 <- function(id, display_cancer_type, display_gene, display_pcx, display_pcy, pca_annotation_na_rm, rv) {
  moduleServer( id,
	function(input, output, session) {
		output$showPlot <- renderPlotly({
			validate(need(display_cancer_type(), message = FALSE), 
				   need(display_gene(), message =  FALSE))

			# use pre-processed data to make PCA plot
			out_plot <- make_pca_plot(
						pca_annotation = display_gene(),
						pcx = display_pcx(),
						pcy = display_pcy(), 
						scree_df = rv[["scree"]][[display_cancer_type()]], 
						pca_df= rv[["pca"]][[display_cancer_type()]],
						gene_full_name = names(rv[["formatted_gene_list"]])[rv[["formatted_gene_list"]] == display_gene()],
						cancer_full_name = names(rv[["input_cancer_type"]])[rv[["input_cancer_type"]] == display_cancer_type()],
						col_label = rv[["exprs_unit"]],
						discrete_cols = discrete_cols,
						remove_na = pca_annotation_na_rm())
			ggplotly(out_plot) 
		})
	}
  )
}


## PCA plot2 (colour by selected annotation)
pcaPlot2 <- function(id, display_cancer_type, display_gene, display_pcx, display_pcy, display_pca_annotation, pca_annotation_na_rm, pca_annotation_figLegend, rv) {
  moduleServer(id,
	function(input, output, session) {
		output$showPlot <- renderPlotly({
			validate(need(display_cancer_type(), message = FALSE), 
				   need(display_gene(), message =  FALSE),
				   need(display_pca_annotation(), message =  FALSE))

			# Organise current annotation
			sampleAnnot_pca <- add_sampleAnnot_pca[add_sampleAnnot_pca %in% colnames(rv[["pca"]][[display_cancer_type()]])]
			
			# use pre-processed data to make PCA plot
			out_plot <- make_pca_plot(
				pca_annotation = display_pca_annotation(),
				pcx = display_pcx(),
				pcy = display_pcy(), 
				scree_df = rv[["scree"]][[display_cancer_type()]], 
				pca_df = rv[["pca"]][[display_cancer_type()]],
				gene_full_name = names(sampleAnnot_pca)[sampleAnnot_pca == display_pca_annotation()],
				cancer_full_name = names(rv[["input_cancer_type"]])[rv[["input_cancer_type"]] ==  display_cancer_type()],
				col_label = "",
				discrete_cols = discrete_cols,
				remove_na = pca_annotation_na_rm()) 
			# toggle figure legend
			if(pca_annotation_figLegend()){
				ggplotly(out_plot) %>% layout(showlegend = F)
			} else {
				ggplotly(out_plot)
			}		})
	}
  )
}

## Additional: Scree plot
pcaScree <- function(id, display_cancer_type, display_gene, rv) {
  moduleServer(
	id,
	function(input, output, session) {
		output$showPlot <- renderPlotly({
			validate(need(display_cancer_type(), message = FALSE), 
				   need(display_gene(), message =  FALSE))

			# use pre-processed data to make PCA plot
			out_plot <- make_scree_plot(
						gene = display_gene(),
						scree_df = rv[["scree"]][[display_cancer_type()]], 
						gene_full_name = names(rv[["formatted_gene_list"]])[rv[["formatted_gene_list"]] == display_gene()],
						cancer_full_name = names(rv[["input_cancer_type"]])[rv[["input_cancer_type"]] == display_cancer_type()])
			ggplotly(out_plot) 
		})
	}
  )
}

## Additional: Gene expression and sample annotation_plot
pcaGeneAnnotation <- function(id, display_cancer_type, display_gene, display_pca_annotation, display_pcy, pca_annotation_na_rm, rv) {
  moduleServer(id,
	function(input, output, session) {
		output$showPlot <- renderPlotly({
			validate(need(display_cancer_type(), message = FALSE), 
				   need(display_gene(), message =  FALSE),
				   need(display_pca_annotation(), message =  FALSE))

			# Organise current annotation
			sampleAnnot_pca <- add_sampleAnnot_pca[add_sampleAnnot_pca %in% colnames(rv[["pca"]][[display_cancer_type()]])]
			
			# use pre-processed data to make PCA plot
			out_plot <- make_gene_annotation_plot(
				pca_annotation = display_pca_annotation(),
				scree_df = rv[["scree"]][[display_cancer_type()]], 
				pca_df = rv[["pca"]][[display_cancer_type()]],
				annotation_full_name = names(add_sampleAnnot_pca)[add_sampleAnnot_pca == display_pca_annotation()],
				input_display_gene = display_gene(),
				input_display_cancer_type = display_cancer_type(),
				gene_full_name = names(rv[["formatted_gene_list"]])[rv[["formatted_gene_list"]] == display_gene()],
				cancer_full_name = names(rv[["input_cancer_type"]])[rv[["input_cancer_type"]] ==  display_cancer_type()],
				input_display_pcy = display_pcy(),
				#col_label = "",
				#discrete_cols = discrete_cols,
				remove_na = pca_annotation_na_rm()) 

			ggplotly(out_plot)
		})
	}
  )
}


## Top gene contribution on each PC
pcaTopGenes <- function(id, display_cancer_type, pca_annotation_figLegend, rv) {
  moduleServer(
	id,
	function(input, output, session) {
		output$showPlot <- renderPlotly({
			validate(need(display_cancer_type(), message = FALSE))

			# use pre-processed data to make PCA plot
			out_plot <- make_topGeneContribution_plot(
						scree_df = rv[["scree"]][[display_cancer_type()]], 
						contribution_df = rv[["top_contributions"]][[display_cancer_type()]], 
						top_gene_contribution = top_gene_contribution_num)
			ggplotly(out_plot)  %>%
			     plotly::layout(legend=list(x=0, 
                                 xanchor='left',
                                 yanchor='bottom',
                                 orientation='h'))  %>%
				layout(xaxis = list(side ="top" ) )	
		})
	}
  )
}


		
