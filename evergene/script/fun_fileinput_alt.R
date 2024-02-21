# based on the Shiny fileInput function
fileInput_alt <- function(inputId, label = NULL, labelIcon = NULL, multiple = FALSE, 
                       accept = NULL, width = NULL, progress = TRUE, ...) {
  # add class fileinput_2 defined in UI to hide the inputTag
  inputTag <- tags$input(id = inputId, name = inputId, type = "file", 
                         class = "fileinput_2")
  if (multiple) 
      inputTag$attribs$multiple <- "multiple"
  if (length(accept) > 0) 
      inputTag$attribs$accept <- paste(accept, collapse = ",")

  div(..., style = if (!is.null(width)) paste0("width: ", validateCssUnit(width), ";"), 
    inputTag,
    # Customise action button
    tags$label(`for` = inputId, div(icon(labelIcon), label, 
                class = "btn btn-default action-button",
				style="padding-top:3px; padding-bottom:3px; padding-left:3px; padding-right:3px;")),
    # Progress bar option
    if(progress)
      tags$div(id = paste(inputId, "_progress", sep = ""), 
        class = "progress shiny-file-input-progress", 
        tags$div(class = "progress-bar")
      )
  )
}     