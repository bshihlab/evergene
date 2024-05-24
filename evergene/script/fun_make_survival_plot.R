##create function START
make_survival_plot <- function(cancer_full_name, 
								gene, 
								gene_full_name,
								labelled_data, 
								plot_cols = c("#89CFF0", "#cc210a"), 
								threshBottom=33.3,
								threshTop=66.7,
								event = "PFI"){ 

	# Organise output object
	out_list <- list(plot_data = NULL, Cox_stats = NULL, Log_rank_stats = NULL, surv_plot = NULL)
    stats_result_vec <- c(
        "gene" = gene,
        "gene_full_name" = NULL,
		"cancer" = cancer_full_name,
        "p.value" = NULL,
        "totalDeath" = NULL,
        "numDeath.highExprs" = NULL,
        "numDeath.lowExprs" = NULL,
        "meanTpm.highExprs" = NULL,
        "meanTpm.lowExprs" = NULL)
	
	
	# remove samples without survival data
	labelled_data <- labelled_data[!is.na(labelled_data[,event]),]
	labelled_data <- labelled_data[!is.na(labelled_data[, paste0(event, ".time")]),]

	# make a current data frame for plotting
	surv_df <- survival_grouping(
			barcode = labelled_data$barcode,
			exprs = labelled_data[[gene]],
			#exprs = log2(labelled_data[[gene]]+0.1),			
			survival_status = as.numeric(labelled_data[,event]), 
			survival_time = as.numeric(labelled_data[,paste0(event, ".time")]),
			threshBottom = threshBottom , 
			threshTop = threshTop)
	# Turn the expression value into log for cox regression if the minimum value is >0
	if(min(surv_df$exprs)>=0){
		surv_df$exprs <- log(surv_df$exprs + 0.1)
	} 


	surv_df_2cat <- surv_df[surv_df$exprs_grp %in% c("Low.exprs", "High.exprs"), ]
	#surv_df$gene <- gene
	#surv_df$gene_full_name <- gene_full_name
	
	legendHigh <- "High.exprs"
	legendLow  <- "Low.exprs"

	# fill in statistics in the output vector
	stats_result_vec["totalDeath"] <- sum(surv_df_2cat$survival_status, na.rm=TRUE)
	stats_result_vec["numDeath.highExprs"] <- sum(surv_df_2cat$survival_status[surv_df_2cat$exprs_grp %in% "High.exprs"])
	stats_result_vec["numDeath.lowExprs"] <- sum(surv_df_2cat$survival_status[surv_df_2cat$exprs_grp %in% "Low.exprs"])
	stats_result_vec["meanTpm.highExprs"] <- mean(surv_df_2cat$exprs[surv_df_2cat$exprs_grp %in% "High.exprs"])
	stats_result_vec["meanTpm.lowExprs"] <-  mean(surv_df_2cat$exprs[surv_df_2cat$exprs_grp %in% "Low.exprs"])

	## Survival fitted as continuous
	# code adapted from http://www.sthda.com/english/wiki/cox-proportional-hazards-model
	# http://www.sthda.com/english/wiki/cox-model-assumptions
	surv_cox_results <- tryCatch({
		#surv_df_cox <- surv_df[!is.na(surv_df$exprs),]
		res.cox <- survival::coxph(survival::Surv(survival_time, survival_status) ~ exprs, data = surv_df)
		x <- summary(res.cox)
		p.value<-signif(x$wald["pvalue"], digits=3)
		wald.test<-signif(x$wald["test"], digits=3)
		beta<-signif(x$coef[1], digits=3);	#coeficient beta
		HR <-signif(x$coef[2], digits=3);	#exp(beta)
		HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
		HR.confint.upper <- signif(x$conf.int[,"upper .95"], 3)
		HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
		res<-c(gene, cancer_full_name, beta, HR, wald.test, p.value)
		names(res) <- c("gene", "cancer", "beta", "HR (95% CI for HR)", "wald.test", "p.value")
		res
	}, error = function(e) {
		res <- c(NA, NA, NA, NA, NA, NA, NA)
		names(res) <-  c("gene",  "cancer", "beta", "HR (95% CI for HR)", "wald.test", "p.value")
		return(res)
	})

	## calculate 2 catogery surivial
	surv_results <- tryCatch({
		tabSurv <- survival::survdiff(data = surv_df_2cat, survival::Surv(survival_time, survival_status)  ~ exprs_grp)
		D1 <- tabSurv$obs[1]
		D2 <- tabSurv$obs[2]
		E1 <- tabSurv$exp[1]
		E2 <- tabSurv$exp[2]
		HR <- (E1/E2)/(D1/D2)
		SE_lnHR = sqrt(1/E1 + 1/E2)
#		HR <-(tabSurv$obs[2]/tabSurv$exp[2])/(tabSurv$obs[1]/tabSurv$exp[1])
	
		L = log(HR)
		HR.confint.lower <- signif(exp(L - 1.96*SE_lnHR), 3)
		HR.confint.upper <- signif(exp(L + 1.96*SE_lnHR), 3)
		tabSurv_chis <- unlist(tabSurv)$chisq
		pval <- as.numeric(1 - pchisq(abs(tabSurv$chisq), df = 1))
		HR <- paste0(signif(HR, digits=3), " (", HR.confint.lower, "-", HR.confint.upper, ")")
		res<-c(pval, HR)
		names(res) <- c("logrank_p", "logrank_hr")
		res
	}, error = function(e) {
		return(NA)
	})
	
	stats_result_vec["p.value"] <-  surv_results["logrank_p"]
		

	## make plots
	# open new plot

	fit <- survival::survfit(data = surv_df_2cat, survival::Surv(survival_time, survival_status)  ~ exprs_grp)
	
	#plot.new() 
	#dev.control("enable")	
#	par(mar = c(2, 2, 9, 2))
	
	# make plot
	pretty_p_val <- prettyNum(surv_results["logrank_p"])
	titlePlot <-
		paste0(
			  cancer_full_name, "\n",
			  gene_full_name, "\n",
			  "Kaplan-Meier Survival curve\n",
			  "Log-rank test [ p.val=", pretty_p_val, ",\n    HR(95%CI)=" , surv_results["logrank_hr"], " ]\n" ,
			  "Cox regression [ p.val=", prettyNum(surv_cox_results["p.value"]), ", beta=", surv_cox_results["beta"],",\n    HR(95%CI)=" , surv_cox_results["HR (95% CI for HR)"], " ]")
#	plot(
#	  fit,
#		col = plot_cols,
#		xlab = "Days",
#		ylab = "Progression-free probability",
#		lwd = 4, main=titlePlot
#	)
	
#	#title(titlePlot, adj=0, font.main = 1)

#	legend(
#		"topright", bty = "n",
#		legend = c(legendLow, legendHigh),
#		col = plot_cols,
#		text.col = "#000000",
#		pch = 15 )
	
	#prnt(paste0("Loop: Surivial analysis : gene loop: ", gene, ":save plots"))
	# save the recorded plot
	#out_list[["surv_plot"]] <- recordPlot()
	
	theme <- theme_bw() + theme(panel.grid=element_blank(),legend.background = element_blank()) 
		   
		   
	p <- survminer::ggsurvplot(fit = fit,
							data = surv_df_2cat,
							#risk.table = TRUE,
							legend =  c(0.8,0.9),
							legend.title = "",
							legend.labs = c("Low.exprs", "High.exprs"),
							palette = c("#89CFF0", "#cc210a"), 
							surv.scale="percent",
							ggtheme = theme
							)
	p <- p + ggtitle(titlePlot) + xlab(paste0(clinical_event[event], " duration (day)")) + ylab(paste0("Probability of ", clinical_event[event]))  

	out_list[["surv_plot"]] <- ggpar(p, 
								  font.main = c(14),
								  font.x = c(14),
								  font.y = c(14),
								  font.caption = c(14), 
								  font.legend = c(14), 
								  font.tickslab = c(14))
	#dev.off()
	


	#### reorganise output
    # Add cancer name to the table
	out_list[["Log_rank_stats"]] <- stats_result_vec
	out_list[["Cox_stats"]] <- surv_cox_results
	out_list[["plot_data"]] <- surv_df
	
	
	return(out_list)
  
}

