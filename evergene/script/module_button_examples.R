#### Example buttons
exampleButton <- function(id, label) {
	ns <- NS(id)
	tagList(
		actionButton(ns("giveExample"), label, 
			style="color: #000000; background-color: #ffffff; border-color: #DFD7CA; padding:6px") 
	)
}

exampleButtonServer <- function(id, selectedGenes, selectedCancer, parent_session) {
  moduleServer(
	id,
	function(input, output, session) {
	  observeEvent(input$giveExample, {
		updateTextAreaInput(session = parent_session, 'gene_string', "", value=selectedGenes)
		updatePickerInput(session = parent_session, "cancer_type_list", "Select Cancer Type", options = list('actions-box' = TRUE), selected=selectedCancer , list_of_cancer_types)
	  })
	}
  )
}
