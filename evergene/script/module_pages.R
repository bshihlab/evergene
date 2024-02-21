### each of the main pages
# Analysis
pages_analysis <- tabPanel("Analysis", 
								sidebarLayout(
									shinyUI_Panel_inupt, 
									shinyUI_Panel_main
								) # Sidebar layout
                           )
# About
pages_about <- tabPanel("About",
		fluidRow(
			column(width = 4, wellPanel(
			HTML("Created by A. Kennedy, E. Richardson and B. Shih. <br/> <br/> 
					This tool was created with funding from North West Cancer Research."), 
		  ), # Well panel
		  column(width = 6, br(), 
				 img(src='nwcr.png', height="100%", width="100%")), 
		  column(width = 6, br(), 
				 img(src='lu.png', height="100%", width="100%"))
		  ), # Column
		  
		  column(width = 8, wellPanel(
		  #  HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/Ka2pWqXS1WA" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>')
		  ) # Well panel
		  ) # Column
		)  # Fluid row
	) # Tab panel

# Contact us
pages_contactUs <- tabPanel("Contact us")

# Source page
pages_source <- tabPanel("Source code", 
	HTML('<script src="https://emgithub.com/embed-v2.js?target=https%3A%2F%2Fgithub.com%2Fellaintheclouds%2FNWCR%2Fblob%2Fmain%2FPCA_function.R&style=default&type=code&showBorder=on&showLineNumbers=on&showFileMeta=on&showFullPath=on&showCopy=on"></script>'))