#setwd("C:/Users/shihb/OneDrive - Lancaster University/work/project/2023/20230816_nwcr/evergene")

# Set up-------------------------------------------------------------------------
source("script/load_libraries.R")
source("script/load_variables.R")
source("script/fun_input_format.R")
source("script/fun_show_invalid_genes.R")
source("script/fun_fileInput_alt.R")
source("script/fun_make_pca_plot.R")
source("script/fun_make_pca_plot_multiGene.R")
source("script/fun_make_gene_project_count_plot.R")
source("script/fun_make_scree_plot.R")
source("script/fun_make_geneContribution_plot.R")
source("script/fun_make_survival_plot.R")
source("script/fun_make_topGeneContribution_plot.R")
source("script/fun_make_gene_annotation_plot.R")
source("script/fun_make_corr_plot.R")
source("script/fun_survival_grouping.R")
source("script/fun_save_pca.R")
source("script/fun_save_survival.R")
source("script/fun_userInput_load.R")
source("script/fun_userInput_pca_topGenes.R")
source("script/module_button_examples.R")
source("script/module_button_analyse.R")
source("script/module_pca.R")
source("script/module_pca_download.R")
source("script/module_survival.R")
source("script/module_survival_download.R")
source("script/module_correlation.R")
source("script/module_correlation_download.R")
source("script/module_shinyUI_style.R")
source("script/module_shinyUI_components.R")
source("script/module_shinyUI.R")
source("script/module_pages.R")

# UI----------------------------------------------------------------------------
ui <- fluidPage(theme = shinytheme("sandstone"),         
	shinyUI_taglist,
	navbarPage(
			   title=div(img(src="evergene_grey.png", height="50", aOlign="left"), "",  class = "dropdown"),
			   pages_analysis, 
			   pages_about,
			   #pages_contactUs, 
			   #pages_source
	) # Navbar page
) # Fluid page


# Server------------------------------------------------------------------------
server <- function(input, output, session){
	
	#### Before analyse -------------------------------------------------------
	## Observer helpers to populate hints
	observe_helpers()
	
	
	
	### Import user-uploaded input data
	# Show buttons
#	renderTable(input$upload_count)
#	renderTable(input$upload_annotation)
	output$upload_log_transform <- renderText({input$upload_log_transform})
	output$upload_raw_count <- renderText({input$upload_raw_count})

	
	## Auto fill buttons for example cases
	exampleButtonServer("example_case1", selectedGenes = "TP63\nSOX2\nGATA4\nGATA6\nPHB", selectedCancer = c("TCGA-ESCA"), parent_session = session)
	exampleButtonServer("example_case2", selectedGenes = "ESR1\nPGR\nERBB2\nCENPN\nBRCA1\nBRCA2", selectedCancer = c("TCGA-BRCA"), parent_session = session)


	## Initiate organised data values
	organised_data <- reactiveValues(
									"dataSource" = "TCGA",
									"showPanel" = FALSE, 
									"exprs_unit" = "TPM", 
									"computingInProgress" = FALSE,
									"input_cancer_type" = NULL, 
									"formatted_gene_list" = NULL, 
									"tpm" = list(), 
									"pca" = list(),
									"scree" = list(),
									"survival" = list(),
									"top_contributions" = list(),
									"max_pc" = num_pcs,
									"dropdown_pca" = list(),
									"dropdown_survival" = list(),
									"dropdown_correlation" = list(),
									"display_message" = "")

	### Output message to show the update on which project is selectd
#	input_project_update <- reactiveVal()
#	input_project_update("Please select at least 1 cancer project or upload data & annotation")
	
	
	## Display barplot/message showing genes found/selected, genes not found, and projects
	# Make reactive plots
#	gene_project_count <- reactive({
#		validate(need(input$cancer_type_list, message = "Please enter cancer types."), need(input$gene_string, message = "Please enter genes.")) # Validate needs
#		make_gene_project_count_plot(input$gene_string, input$cancer_type_list)
#	})
#
#	# display gene-cancer combination count plots
#	output$display_pca_plot1 <- renderPlotly({gene_project_count()})
#	output$display_survival_plot1 <- renderPlotly({gene_project_count()})


	# Show invalid genes
	output$invalid_gene_message_display <- renderText({
		validate(need(input$gene_string, message = FALSE)) # Validate needs
		# show non matching genes if the current dataSourece is TCGA
		if(organised_data[["dataSource"]] == "TCGA"){
			paste(show_invalid_genes(input$gene_string, input$cancer_type_list))
		}
	})
#	output$input_project_update <- renderText({
#		input_project_update()
#	})

	###################################################################
	#### Events that occur after pressing "Analyse"  -------------------------------------------------------

#print("step1")
	## IMPORT DATA
	# After clicking on "Analyse", this module takes in the current input and updates the reactive values "organise_data"
	analyseButtonServer("button_compute", 
		input_cancer_type = reactive({input$cancer_type_list}), 
		input_gene_string = reactive({input$gene_string}), 
		user_input_show = reactive({ input$user_input_show }),
		upload_count = reactive({input$upload_count}), 
		upload_annotation = reactive({input$upload_annotation}),
		upload_gene = reactive({input$upload_gene}),
		upload_log_transform = reactive({input$upload_log_transform}),
		upload_raw_count = reactive({input$upload_raw_count}),
		rv = organised_data
#		rv_message = input_project_update
		) 
		

	#### Show main panels after analysis
	#output$computed <- renderText(organised_data$input_cancer_type)
	output$panelStatus <- reactive({ organised_data$showPanel })
	outputOptions(output, "panelStatus", suspendWhenHidden = FALSE)

	### Turn on/off user input
	# turn on/off the user input boxes
	output$user_input_show <- reactive({input$user_input_show})
	outputOptions(output, "user_input_show", suspendWhenHidden = FALSE)
	

#print("step2")

	#### Main1 : dropdown menu 
	# dropdown menu (selected cancer, gene, pc_x, pc_y)
	# Display the dropdown accordingly
	output$display_cancer_type <- renderUI({selectInput("display_cancer_type", "Cancer type", organised_data$input_cancer_type, selected = organised_data$input_cancer_type[1])})
	output$display_gene <- renderUI({selectInput("display_gene", "Gene", organised_data$formatted_gene_list, selected = organised_data$formatted_gene_list[1])})		
	output$display_pcx <- renderUI({selectInput("display_pcx", "PC (x-axis)", c(1:organised_data$max_pc), selected = 1)})
	output$display_pcy <- renderUI({selectInput("display_pcy", "PC (y-axis)", c(1:organised_data$max_pc), selected = 2)})
	
	#### Main2: PCA gene contribut (top panel)
	pcaGeneContributionServer("display_gene_contribution", display_cancer_type= reactive({input$display_cancer_type}), display_gene=reactive({input$display_gene}), rv = organised_data)

	#### Main3: dropdown menu 
	listenChange <- reactive({
		list(input$display_cancer_type,
			organised_data$computingInProgress)
	})
	observeEvent(listenChange(), {
		validate(need(input$display_cancer_type, message = "Please use the side panel to enter input cancer types."), 
				need(input$display_gene, message = "Please use the side panel to enter input genes.")) # Validate needs

		# PCA dropdown : sample annotation
		output$display_pca_annotation <- renderUI({selectInput("display_pca_annotation", "Sample annotation", organised_data[["dropdown_pca"]][[input$display_cancer_type]], selected = organised_data[["dropdown_pca"]][[input$display_cancer_type]][1])})

		# Correlation dropdown
		output$display_correlation_annotation1 <- renderUI({selectInput("display_correlation_annotation1", "X-axis", organised_data[["dropdown_correlation"]][[input$display_cancer_type]], selected = organised_data[["dropdown_correlation"]][[input$display_cancer_type]][1])})
		output$display_correlation_annotation2 <- renderUI({pickerInput("display_correlation_annotation2", "Y-axis (accepts multiple selection)", options = list('actions-box' = TRUE), organised_data[["dropdown_correlation"]][[input$display_cancer_type]], multiple = TRUE, selected = organised_data[["dropdown_correlation"]][[input$display_cancer_type]][2])})

		# Survival dropdown : events
		output$display_survival_event <- renderUI({selectInput("display_survival_event", "Events", clinical_events_all[clinical_events_all %in% colnames(organised_data[["survival"]][[input$display_cancer_type]])] , selected = "PFI")})

		#  Survival dropdown : sample annotation
		output$display_survival_annotation <- renderUI({selectInput("display_survival_annotation", "Sample annotation", organised_data[["dropdown_survival"]][[input$display_cancer_type]], selected = organised_data[["dropdown_survival"]][[input$display_cancer_type]][1])})
	})
#print("step3")
	
	#### Main3: PCA (plots)
	# Plot1: coloured by gene expression
	pcaPlot1("display_pca_plot1", 
		display_cancer_type = reactive({input$display_cancer_type}), 
		display_gene = reactive({input$display_gene}), 
		display_pcx = reactive({input$display_pcx}), 
		display_pcy = reactive({input$display_pcy}), 
		pca_annotation_na_rm = reactive({input$pca_annotation_na_rm}), 
		rv = organised_data)
		
	# Plot2: coloured by gene expression
	pcaPlot2("display_pca_plot2", 
		display_cancer_type = reactive({input$display_cancer_type}), 
		display_gene = reactive({input$display_gene}), 
		display_pcx = reactive({input$display_pcx}), 
		display_pcy = reactive({input$display_pcy}), 
		display_pca_annotation = reactive({input$display_pca_annotation}), 
		pca_annotation_na_rm = reactive({input$pca_annotation_na_rm}), 
		pca_annotation_figLegend = reactive({input$pca_annotation_figLegend}),
		rv = organised_data)
		
	# Additional : Screeplot
	pcaScree("display_scree_plot", 
		display_cancer_type = reactive({input$display_cancer_type}), 
		display_gene = reactive({input$display_gene}), 
		rv = organised_data)

	# Additional : Gene exp and Sample annotation
	pcaGeneAnnotation("display_geneAnnotation_plot", 
		display_cancer_type = reactive({input$display_cancer_type}), 
		display_gene = reactive({input$display_gene}), 
		display_pca_annotation = reactive({input$display_pca_annotation}), 
		display_pcy = reactive({input$display_pcy}), 
		pca_annotation_na_rm = reactive({input$pca_annotation_na_rm}), 
		rv = organised_data)
		
	#  Top gene contribution on each PC
	pcaTopGenes("display_tile_plot", 
		display_cancer_type = reactive({input$display_cancer_type}),
		pca_annotation_figLegend = reactive({input$pca_annotation_figLegend}),		
		rv = organised_data)
#print("step4")
	
	#### Main3: Survival (plots)
	# Plot1: coloured by gene expression
	survivalPlot1("display_survival_plot1", 
		display_cancer_type = reactive({input$display_cancer_type}), 
		display_gene = reactive({input$display_gene}), 
		survival_percent_lower = reactive({input$survival_percent_lower}), 
		survival_percent_upper = reactive({input$survival_percent_upper}), 
		display_survival_event = reactive({input$display_survival_event}), 
		rv = organised_data)

	# Plot2: coloured by annotation/gene expression
	survivalPlot2("display_survival_plot2", 
		display_cancer_type = reactive({input$display_cancer_type}), 
		display_gene = reactive({input$display_gene}), 
		survival_percent_lower = reactive({input$survival_percent_lower}), 
		survival_percent_upper = reactive({input$survival_percent_upper}), 
		display_survival_event = reactive({input$display_survival_event}),
		display_survival_annotation = reactive({input$display_survival_annotation}),
		display_pcy = reactive({input$display_pcy}),
		survival_annotation_na_rm = reactive({input$survival_annotation_na_rm}),
		survival_annotation_figLegend = reactive({input$survival_annotation_figLegend}),
		rv = organised_data)
#print("step5")

	#### Main3: Correlation (plots)
	# Plot1: coloured by gene expression
	correlationPlot1("display_corr_plot", 
		display_cancer_type = reactive({input$display_cancer_type}), 
		display_correlation_annotation1 = reactive({input$display_correlation_annotation1}), 
		display_correlation_annotation2 = reactive({input$display_correlation_annotation2}), 
		correlation_log_gene_exp = reactive({input$correlation_log_gene_exp}), 
		rv = organised_data)
#print("step6")

	#### Main3: Download buttons 
	## PCA
	pcaDownloadPdf("download_pdf_pca_button", 
		rv = organised_data)		
	pcaDownloadZip("download_zip_pca_button", 
		display_pcx = reactive({input$display_pcx}),
		display_pcy = reactive({input$display_pcy}),
		display_pca_annotation = reactive({input$display_pca_annotation}),
		rv = organised_data)
	## Survival
	survivalDownloadPdf("download_pdf_survival_button",
		survival_percent_lower = reactive({input$survival_percent_lower}), 
		survival_percent_upper = reactive({input$survival_percent_upper}), 
		display_survival_event = reactive({input$display_survival_event}),
		rv = organised_data)
	survivalDownloadZip("download_zip_survival_button",
		survival_percent_lower = reactive({input$survival_percent_lower}), 
		survival_percent_upper = reactive({input$survival_percent_upper}), 
		display_survival_event = reactive({input$display_survival_event}),
		rv = organised_data)
	## Correlation
	correlationDownloadPdf("download_pdf_correlation_button", 
		display_correlation_annotation1 = reactive({input$display_correlation_annotation1}), 
		correlation_log_gene_exp = reactive({input$correlation_log_gene_exp}),
		rv = organised_data)
	correlationDownloadZip("download_zip_correlation_button", 
		display_correlation_annotation1 = reactive({input$display_correlation_annotation1}), 
		correlation_log_gene_exp = reactive({input$correlation_log_gene_exp}),
		rv = organised_data)

} # Server


# Run---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
