userInput_pca_topGenes <- function(current_df, gene_annotation){
	row.names(current_df) <- current_df$PC
	current_df$PC=NULL
	current_df$contribution=NULL
	current_df <- t(current_df)
	current_df_abs <- abs(current_df)
	current_df_rank <- as.data.frame(apply(current_df_abs, 2, order, decreasing=TRUE))
	# get table for top_genes
	# loop through each column
	top_genes <- 100
	top_gene_contributions <- list()
	top_gene_ids <- list()
	for(idx in 1:ncol(current_df_rank)){
	  current_colname <- colnames(current_df_rank)[idx]
	  top_gene_contributions[[current_colname]] <- current_df[,idx][current_df_rank[,idx]]
	  top_gene_contributions[[current_colname]] <- top_gene_contributions[[current_colname]][1:top_genes]
	  names(top_gene_contributions[[current_colname]]) <- 1:top_genes
	  top_gene_ids[[current_colname]] <- row.names(current_df)[current_df_rank[,idx]]
	  top_gene_ids[[current_colname]] <- top_gene_ids[[current_colname]][1:top_genes]
	}

	top_gene_contributions <- do.call(cbind, top_gene_contributions)
	top_gene_ids <- do.call(cbind, top_gene_ids)
	
	# get readable gene names
	top_gene_names <- apply(top_gene_ids, 2, function(x){matched_val <- gene_annotation[gene_annotation$gene_id %in% x,]; matched_val$gene_name[match(x, matched_val$gene_id)]})
	
	# reorder columns
	top_gene_names <- as.data.frame(top_gene_names[,paste0("PC", 1:ncol(current_df_rank))])
	top_gene_contributions <- as.data.frame(top_gene_contributions[,paste0("PC", 1:ncol(current_df_rank))])
	
	# wide to long
	organised_val <- list()
	for(i in 1:ncol(current_df_rank)){
		current_PC <- paste0("PC", i)
		organised_val[[current_PC]] <- data.frame(PC = current_PC, gene_rank = 1:top_genes, gene = top_gene_names[[current_PC]], contribution = top_gene_contributions[[current_PC]])
	}
	organised_val <- do.call(rbind,organised_val)
	return(organised_val)
}