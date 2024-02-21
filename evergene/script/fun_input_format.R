# 1 Function--------------------------------------------------------------------
input_format <- function(gene_list, cancer_list, in_gene_annotation = tcga_gene_annotation, combination_restriction = pca_plot_combination_restriction ){
	# split by comma, space or new line
	gene_list <- strsplit(gene_list, split = "\n| |,|;")[[1]] # Split the user input into a list
	gene_list <- gene_list[!gene_list == ""] # Remove "" from gene vector (introduced when user inputs multiple new lines, spaces or comma)
	gene_list <- unique(gene_list)
	
	# initiate lists
	formatting_output <- list("rejected_gene_list" = vector(), "formatted_gene_list" =list())
	idx <- 1
	for (g in gene_list){
	g_upper <- toupper(g)

	in_gene_annotation$gene_name_upper <- toupper(in_gene_annotation$gene_name)

	if (g_upper %in% in_gene_annotation$gene_id | g_upper %in% in_gene_annotation$gene_id_stable |g_upper %in%  in_gene_annotation$gene_name_upper){
	  
	  gene_annotation_relavent <- in_gene_annotation[in_gene_annotation$gene_id == g_upper | in_gene_annotation$gene_id_stable == g_upper | in_gene_annotation$gene_name == g_upper | in_gene_annotation$gene_name_upper == g_upper, c("gene_id", "gene_name")]
	  # create a merged_name from gene_name and gene_id if the gene_name value is not ""
	  gene_annotation_relavent$merged_name <- ifelse((gene_annotation_relavent$gene_name == "" | gene_annotation_relavent$gene_name == gene_annotation_relavent$gene_id), 
														gene_annotation_relavent$gene_id, 
													   paste0(gene_annotation_relavent$gene_name, ": ", gene_annotation_relavent$gene_id))
	  formatting_output[["formatted_gene_list"]][[idx]] <-   gene_annotation_relavent
	  
	  idx <- idx+1
	} # If

	else {
	  formatting_output[["rejected_gene_list"]] <- c(formatting_output[["rejected_gene_list"]], g)
	} # Else
	} # For

	# Make dataframe from list
	genes_df <-  do.call(rbind, formatting_output[["formatted_gene_list"]])
	genes_df <- unique(genes_df) 	# Keep only unique entries
	
	# check the combination restriction. Return gene list shorter than the input gene list so that a maximum of 500 gene- cancer combinations are inputted
	max_gene_keep <- floor(combination_restriction / length(cancer_list))
	
	# reduce the gene list to include only gene numbers that will keep under the thresholds (i.e. if 25 cancer types were selected, only the first 20 genes will be included)
	if(length(genes_df) > max_gene_keep){
		genes_df <- genes_df[1:max_gene_keep,]
	}
	
	 
	formatting_output[["formatted_gene_list"]] <- genes_df

		
	# Add another element in the list so it has a named-vector with 
	# gene IDs [names as merged_name] for displaying in shinyapp drop down
	formatting_output[["formatted_gene_list_named"]] <- formatting_output[["formatted_gene_list"]]$gene_id
	names(formatting_output[["formatted_gene_list_named"]]) <- formatting_output[["formatted_gene_list"]]$merged_name
	formatting_output[["formatted_gene_list_named"]] <- formatting_output[["formatted_gene_list_named"]][!is.na(formatting_output[["formatted_gene_list_named"]])]

	
	# Return output
	return(formatting_output)

} # function