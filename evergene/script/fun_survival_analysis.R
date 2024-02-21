##create function START
survival_function <- function(cancer_full_name, 
								gene_list, tpm, 
								sample_annotation, 
								gene_annotation, 
								plot_cols = c("#89CFF0", "#cc210a"), 
								threshTop=0.67, 
								threshBottom=0.33,
								event = "PFI"){ 

	# organise data
	labelled_data <- merge(sample_annotation, tpm, by.x="barcode", by.y=0)

	
	# remove duplicated samples (keeping the first one)
	labelled_data <- labelled_data[!duplicated(labelled_data$bcr_patient_barcode),]
	
	# remove samples without survival data
	labelled_data <- labelled_data[!is.na(labelled_data[,event]),]
	labelled_data <- labelled_data[!is.na(labelled_data[, paste0(event, ".time")]),]

	# organise output object
	out_plot_list <- list()
	stats_result_mx <- matrix(0, length(gene_list), 7)
    colnames(stats_result_mx) <- c(
        "gene",
        "pvalue",
        "totalDeath",
        "numDeath.highExprs",
        "numDeath.lowExprs",
        "meanTpm.highExprs",
        "meanTpm.lowExprs")
	
	
	#prnt("# Surivial analysis : begin gene loop")
	#prnt(gene_list)
	# loop through each gene
	for(idx in 1:length(gene_list)){
		gene <- gene_list[idx]
		stats_result_mx[idx, "gene"] <-  gene
		#prnt(paste0("Loop: Surivial analysis : gene loop: ", gene))

		# make a current data frame for plotting
		surv_df <- data.frame(survival_time = as.numeric(labelled_data[,paste0(event, ".time")]), 
							survival_status = as.numeric(labelled_data[,event]), 
							exprs = labelled_data[[gene]],
							barcode = labelled_data$barcode)
				
		# remove na
		surv_df <- surv_df[!is.na(surv_df$survival_status),]
		surv_df <- surv_df[!is.na(surv_df$survival_time),]
		
		# Organise which rows will be kept
		thresholdTop_num <- nrow(surv_df) - round(nrow(surv_df)*threshTop)
		thresholdBottom_num <- round(nrow(surv_df)*threshBottom)
		thresholdNeither_num <- nrow(surv_df) - thresholdTop_num - thresholdBottom_num
		
		# label top and bottom expression groups			
		surv_df <- surv_df[order(surv_df$exprs, surv_df$barcode),]
		surv_df$exprs_grp <- c(rep("Low.exprs", thresholdBottom_num), rep("Mid.exprs", thresholdNeither_num), rep("High.exprs", round(nrow(surv_df)*threshBottom)))

		surv_df_2cat <- surv_df[surv_df$exprs_grp %in% c("Low.exprs", "High.exprs"), ]
		legendHigh <- paste(gene, "High.exprs")
        legendLow  <- paste(gene, "Low.exprs")

		#prnt(paste0("Loop: Surivial analysis : gene loop: ", gene, ":statistics"))
		# fill in statistics
		stats_result_mx[idx, "totalDeath"] <- sum(surv_df_2cat$survival_status, na.rm=TRUE)
		stats_result_mx[idx, "numDeath.highExprs"] <- sum(surv_df_2cat$survival_status[surv_df$exprs_grp %in% "high.exprs"])
		stats_result_mx[idx, "numDeath.lowExprs"] <- sum(surv_df_2cat$survival_status[surv_df$exprs_grp %in% "low.exprs"])
		stats_result_mx[idx, "meanTpm.highExprs"] <- mean(surv_df_2cat$exprs[surv_df$exprs_grp %in% "high.exprs"])
		stats_result_mx[idx, "meanTpm.lowExprs"] <-  mean(surv_df_2cat$exprs[surv_df$exprs_grp %in% "low.exprs"])

		# calculate 2 catogery surivial
		surv_results <- tryCatch({
			tabSurv <- survival::survdiff(data = surv_df_2cat, survival::Surv(survival_time, survival_status)  ~ exprs_grp)
			tabSurv_chis <- unlist(tabSurv)$chisq
			surv_results <- as.numeric(1 - pchisq(abs(tabSurv$chisq), df = 1))
		}, error = function(e) {
			return(NA)
		})
		stats_result_mx[idx, "pvalue"] <-  surv_results
		
		#prnt(paste0("Loop: Surivial analysis : gene loop: ", gene, ":make plots"))
		## make plots
		# open new plot
		plot.new() 
		dev.control("enable")
		
		fit <- survival::survfit(data = surv_df_2cat, survival::Surv(survival_time, survival_status)  ~ exprs_grp)
		
		# open new plot
		plot.new() 
		dev.control("enable")
		
		# make plot
		pretty_p_val <- prettyNum(surv_results)
		titlePlot <-
			paste0("Kaplan-Meier Survival analysis\n",
				   "pvalue=", pretty_p_val, "\n" ,
				  cancer_full_name)
		plot(
		  fit,
			col = plot_cols,
			xlab = "Days",
			ylab = "Survival",
			lwd = 4
		)
		
		title(titlePlot, adj=0, font.main = 1)

		legend(
			"topright", bty = "n",
			legend = c(legendLow, legendHigh),
			col = plot_cols,
			text.col = "#000000",
			pch = 15 )
		
		#prnt(paste0("Loop: Surivial analysis : gene loop: ", gene, ":save plots"))
		# save the recorded plot
		out_plot_list[[gene]] <- recordPlot()
		dev.off()
		
		## Survival fitted as continuous
		# http://www.sthda.com/english/wiki/cox-proportional-hazards-model
		# http://www.sthda.com/english/wiki/cox-model-assumptions
		#res.cox <- survival::coxph(survival::Surv(PFI.time, PFI) ~ exprs, data = surv_df_2cat)
		
	}
	#prnt("# Surivial analysis : end gene loop")

	
	# reorganise output
    stats_result_mx <- as.data.frame(stats_result_mx)

    # Filtering by selected pvalue < 0.01
    rownames(stats_result_mx) <- stats_result_mx$gene
    stats_result_mx[order(stats_result_mx$pvalue, decreasing = FALSE), ]
    
    # Add cancer name to the table
    stats_result_mx$gene_id <- row.names(stats_result_mx)
    stats_result_mx$cancer <- cancer_full_name
    
    # annotate gene name
    stats_result_mx <- merge(gene_annotation[,c("gene_id", "gene_name")], stats_result_mx, by="gene_id")
    stats_result_mx$p.adjust <- p.adjust(stats_result_mx$pvalue, "BH", n=nrow(stats_result_mx)) 
	
	#prnt("# Surivial analysis : end outut reorganisation")
	
	return(list(stats=stats_result_mx, plots=out_plot_list))
  
}

