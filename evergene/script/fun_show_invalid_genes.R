#### Display invalid genes
# Organise plotting dataframe
show_invalid_genes <- function(gene_string, cancer_list){
	formatting_output <- input_format(gene_string, cancer_list)
	genes_notFound <- formatting_output[["rejected_gene_list"]]

	if(length(genes_notFound) > 0){
		return(paste0("Genes not found: ", paste(genes_notFound, collapse=", ")))
	} else {
		return("")
	}
}