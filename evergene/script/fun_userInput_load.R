userInput_load <- function(in_count, in_annotation, in_gene, input_gene_string, in_log_transform, in_raw_count){
	shiny::withProgress(
		message = "Analysing",
		value = 0,
		{

			### valid inputs
			out <- list()
			# check count df
			ext <- tools::file_ext(in_count$name)
			user_count_df <- switch(ext,
			  csv = read.csv(in_count$datapath),
			  tsv = read.delim(in_count$datapath, sep = "\t"),
			  txt = read.delim(in_count$datapath, sep = "\t"),
			  FALSE)
			# check annotation df
			ext <- tools::file_ext(in_annotation$name)
			user_annotation_df <- switch(ext,
			  csv = read.csv(in_annotation$datapath),
			  tsv = read.delim(in_annotation$datapath, sep = "\t"),
			  txt = read.delim(in_annotation$datapath, sep = "\t"),
			  FALSE)
			# check annotation df
			if(isTruthy(in_gene)){
				ext <- tools::file_ext(in_gene$name)
				user_gene_df <- switch(ext,
				  csv = read.csv(in_gene$datapath),
				  tsv = read.delim(in_gene$datapath, sep = "\t"),
				  txt = read.delim(in_gene$datapath, sep = "\t"),
				  FALSE)
			} else {
				user_gene_df <- NULL
			}
			# check if both dataframe exists
			user_count_df_valid <- is.data.frame(user_count_df)
			user_annotation_df_valid <- is.data.frame(user_annotation_df)
			user_gene_df_valid <- is.data.frame(user_gene_df)
			
			# Print progress
			shiny::incProgress(0.25)

			### Organise data if all input present
			if(user_count_df_valid & user_annotation_df_valid ){
				## Organise valid samples
				user_annotation_df <- as.data.frame(user_annotation_df)
				# remove any dashes in column/rownames
				user_annotation_df[,1] <- gsub("-", "_", user_annotation_df[,1])
				rownames(user_annotation_df) <- user_annotation_df[,1]
				colnames(user_count_df) <- gsub("-", "_", colnames(user_count_df))
				
				# remove any duplicated entries for sample annotation
				user_annotation_df <- user_annotation_df[!duplicated(user_annotation_df),]
				
				# check how many sample columns are matched between sample annotation and data
				count_mx <- user_count_df[, colnames(user_count_df)[colnames(user_count_df) %in% rownames(user_annotation_df)]]
				rownames(count_mx) <- user_count_df[,1]
				matched_sample <- ncol(count_mx)
				
				# print error if less than 5 samples
				if(matched_sample < 5){
					out[["message"]] <- "Less than 5 matching samples. Please make sure input file formats are correct. The count matrix and sample annotation have the same sample names, and there are no spaces or special characters (such as -[/!)."
				} else {
					### Organise genes
					# Check if there is gene list input, if not, use the first column as gene_id and gene_name 
					if(user_gene_df_valid){
						col1 <- user_gene_df[,1]
						if(ncol(user_gene_df)>1){col2  <- user_gene_df[,2]} else {col2  <- user_gene_df[,1]}
						user_gene_annotation <- data.frame(gene_id = col1, gene_name = col2, gene_id_stable=col1)
					} else {
						user_gene_annotation <- data.frame(gene_id = user_count_df[,1] , 
													gene_name = user_count_df[,1], 
													gene_id_stable = user_count_df[,1])
					}
					user_gene_annotation <- user_gene_annotation[!duplicated(user_gene_annotation$gene_id),]
					count_mx <- count_mx[row.names(count_mx) %in% user_gene_annotation$gene_id, ]
					matched_gene <- nrow(count_mx)

					# Remove any rows with non-numeric values
					non_numeric_var_idx <- which(apply(count_mx, 1, function(x) sum(is.na(as.numeric(as.character(x)))) > 0))
					if(length(non_numeric_var_idx) > 0){
						count_mx <- count_mx[-non_numeric_var_idx,]
					}
					

					# Remove any rows where there is zero variance
					zero_var_idx <- which(apply(count_mx, 1, function(x) var(x) == 0))
					if(length(zero_var_idx)> 0){
						count_mx <- count_mx[-zero_var_idx,]
					}
					
					# Make sure to convert all values to numeric
					count_mx_row_name <- rownames(count_mx)
					count_mx_col_name <- colnames(count_mx)
					count_mx <- t(apply(count_mx, 1, function(x) as.numeric(as.character(x))))
					count_mx <- as.data.frame(count_mx)
					rownames(count_mx) <- count_mx_row_name
					colnames(count_mx) <- count_mx_col_name
					
					## Calculate from count data
					if(in_raw_count){
						# Normalise data
						x <- DGEList(counts = count_mx)
						keep.exprs <- filterByExpr(x, min.prop=0.2)	# genes expressed in at least 20% of the population
						x <- x[keep.exprs,, keep.lib.sizes=FALSE]
						x <- calcNormFactors(x, method = "TMM")
						count_mx <- cpm(x, log=TRUE)
					} 
					
					
					kept_gene_num <- nrow(count_mx)

					# Note down number of genes and samples are matched
					out[["message"]] <- paste0(matched_sample, " matching samples and ", matched_gene ," genes are found between dataset and annotations. ", kept_gene_num , " with no missing values and > 0 variance are kept.")  


					### Input gene list
					user_gene_formatted <- input_format(input_gene_string, "User input", 
									combination_restriction = pca_plot_combination_restriction, 
									in_gene_annotation = user_gene_annotation)
					formatted_gene_list <- user_gene_formatted[["formatted_gene_list_named"]]
					
					# If input gene not found in the data, use the first gene as placeholder
					if(length(formatted_gene_list)==0){
						first_gene_id <- user_gene_annotation$gene_id[user_gene_annotation$gene_id %in% rownames(count_mx)][1]
						first_gene_name <- user_gene_annotation$gene_name[user_gene_annotation$gene_id %in% rownames(count_mx)][1]
						gene_placeholder <- ifelse(first_gene_name == first_gene_id, first_gene_id, paste0(first_gene_name, ": ", first_gene_id))
						formatted_gene_list <- c(first_gene_id)
						names(formatted_gene_list) <- paste0("Input gene not found: ", first_gene_id, " as place holder")
					}
					
					# add error message about the placeholder
					out[["message"]] <- paste0(out[["message"]] , "\n",  "Gene input not found in the dataset. Using the first gene as a placeholder.")
					
					
					# Print progress
					shiny::incProgress(0.5)					
					
					
					### PCA------------------------------------------------------------------------
					# Process data for PCA

					# the max number of samples (and therefore PCs)
					current_max_pc <- ifelse(matched_sample < num_pcs, matched_sample, num_pcs)		


						
					# Log transform data if requested
					if(in_log_transform){
						# if the lowest value is bigger than zero
						# skip this step if calculated from raw count (log already applied)
						if(in_raw_count){
							count_mx <- count_mx
						} else if(min(count_mx > 0)){
							count_mx <- log10(count_mx)
						} else if (min(count_mx == 0)) {
							count_mx <- log10(count_mx + 1)
						} else {
							out[["message"]] <- paste0(out[["message"]] , "\n",  "Log transformation not performed due to presence of negative values.")
						}
					}
								
					# Run PCA
					set.seed(1) # keep consistent output
					pca_result <- prcomp(t(data.matrix(count_mx)), scale=TRUE) # Carry out PCA on transformed matrix of count data.
					pca_summary <- summary(pca_result)$importance[2,] * 100 # Showing the contribution of each position, to represent its importance.

					contribution_dataframe <- data.frame(PC = paste0("PC", 1:current_max_pc), contribution = pca_summary[1:current_max_pc])

					# Results for Variables (gene contribution)
					res.var <- get_pca_var(pca_result)
					gene_contribution <- res.var$contrib			# Gene contributions to the PCs

					gene_contribution <- res.var$contrib[,1:current_max_pc]
					# add the (+/-)sign for each gene to contribution
					gene_contribution <- as.data.frame(sign(res.var$coord[,1:current_max_pc] )) * as.data.frame(gene_contribution)
					gene_contribution <- t(gene_contribution)

					rownames(gene_contribution) <- paste0("PC", 1:current_max_pc)

					# Print progress
					shiny::incProgress(0.75)					
					
	
					# PC contribution
					scree_df <- merge(contribution_dataframe, gene_contribution, by.x="PC", by.y=0)		# scree PCA
					

					# The x-/y-coordinates to plot for top PCs
					pca_positions <- as.data.frame(pca_result$x)
					pca_positions <- pca_positions[,1:current_max_pc]
					
					pca_data <- as.data.frame(t(count_mx[formatted_gene_list, , drop=FALSE]))
					sample_annotation <- merge(user_annotation_df, pca_data, by.x=colnames(user_annotation_df)[1], by.y=0)
					sample_annotation <- merge(sample_annotation, pca_positions, by.x=colnames(sample_annotation)[1], by.y=0, all.x=TRUE)
					

					out[["pca"]] <- sample_annotation
					out[["max_pc"]] <- matched_sample
					out[["top_contributions"]] <- userInput_pca_topGenes(scree_df, user_gene_annotation[,1:2])
					out[["scree"]] <- scree_df
					out[["formatted_gene_list"]] <- formatted_gene_list
					
					# organise dropdown menu
					sample_annotation_names <- colnames(out[["pca"]])
					names(sample_annotation_names) <- sample_annotation_names
					sample_annotation_names <- sample_annotation_names[2:length(sample_annotation_names)]
					# exclude gene columns for sample annotation 
					sample_annotation_names <- sample_annotation_names[!sample_annotation_names %in% formatted_gene_list]
					out[["dropdown_pca"]] <- sample_annotation_names
					
					# keep only numeric columns for correlation dropdown
					numeric_logic <- sapply(sample_annotation_names, function(x)is.numeric(out[["pca"]][[x]]))
					
					out[["dropdown_correlation"]] <-  c(formatted_gene_list, sample_annotation_names[numeric_logic])
				}
				
			} else {
				missing_message <- vector()
				if (user_count_df_valid == FALSE){ missing_message <- c(missing_message, "count")}
				if (user_annotation_df_valid == FALSE){ missing_message <- c(missing_message, "sample annotation")} 
				out[["message"]] <- paste0("invalid ", paste0(missing_message, collapse=" ,"), " inputs. Please check the help message next to 'User input' switch for details on required inputs or turn off the switch to use TCGA data.")
			}
			
			print("User input imported")
			return(out)
		}
	)
}
