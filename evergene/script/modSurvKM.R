TCGAanalyze_SurvivalKM <- function(
    clinical_patient,
    dataGE,
    Genelist,
    Survresult = FALSE,
    ThreshTop = 0.67,
    ThreshDown = 0.33,
    p.cut = 0.05,
    group1,
    group2,
    cancer_name,
    gene_annotation) {

    #check_package("survival")
    # Check which genes we really have in the matrix
    Genelist <- intersect(rownames(dataGE), Genelist)

    # Split gene expression matrix by the groups
    dataCancer <- dataGE[Genelist, group2, drop = FALSE]
    dataNormal <- dataGE[Genelist, group1, drop = FALSE]
    colnames(dataCancer)  <- substr(colnames(dataCancer), 1, 12)

    cfu <-
        clinical_patient[clinical_patient[, "bcr_patient_barcode"] %in% substr(colnames(dataCancer), 1, 12), ]
    if ("days_to_last_followup" %in% colnames(cfu))
        colnames(cfu)[grep("days_to_last_followup", colnames(cfu))] <-
        "days_to_last_follow_up"
    cfu <-
        as.data.frame(subset(
            cfu,
            select = c(
                "bcr_patient_barcode",
                "days_to_death",
                "days_to_last_follow_up",
                "vital_status"
            )
        ))

    # Set alive death to inf
    if (length(grep("alive", cfu$vital_status, ignore.case = TRUE)) > 0)
        cfu[grep("alive", cfu$vital_status, ignore.case = TRUE), "days_to_death"] <-
        "-Inf"

    # Set dead follow up to inf
    if (length(grep("dead", cfu$vital_status, ignore.case = TRUE)) > 0)
        cfu[grep("dead", cfu$vital_status, ignore.case = TRUE), "days_to_last_follow_up"] <-
        "-Inf"

    cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
    cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]

    followUpLevel <- FALSE

    #FC_FDR_table_mRNA
    tabSurv_Matrix <-
        matrix(0, nrow(as.matrix(rownames(dataNormal))), 8)
    colnames(tabSurv_Matrix) <- c(
        "mRNA",
        "pvalue",
        "Cancer Deaths",
        "Cancer Deaths with Top",
        "Cancer Deaths with Down",
        "Mean Tumor Top",
        "Mean Tumor Down",
        "Mean Normal"
    )

    tabSurv_Matrix <- as.data.frame(tabSurv_Matrix)

    cfu$days_to_death <- as.numeric(as.character(cfu$days_to_death))
    cfu$days_to_last_follow_up <-
        as.numeric(as.character(cfu$days_to_last_follow_up))
    rownames(cfu) <- cfu[, "bcr_patient_barcode"] #mod1

    cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
    cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]

    cfu_complete <- cfu
    ngenes <- nrow(as.matrix(rownames(dataNormal)))
    out_plot_list <- list()
    # Evaluate each gene
    for (i in 1:nrow(as.matrix(rownames(dataNormal))))  {
        cat(paste0((ngenes - i), "."))
        mRNAselected <- as.matrix(rownames(dataNormal))[i]
        mRNAselected_values <-
            dataCancer[rownames(dataCancer) == mRNAselected, ]
        mRNAselected_values_normal <-
            dataNormal[rownames(dataNormal) == mRNAselected, ]
        if (all(mRNAselected_values == 0))
            next # All genes are 0
        tabSurv_Matrix[i, "mRNA"] <- mRNAselected


        # Get Thresh values for cancer expression
        mRNAselected_values_ordered <-
            sort(mRNAselected_values, decreasing = TRUE)
        mRNAselected_values_ordered_top <-
            as.numeric(quantile(as.numeric(mRNAselected_values_ordered), ThreshTop)[1])
        mRNAselected_values_ordered_down <-
            as.numeric(quantile(as.numeric(mRNAselected_values_ordered), ThreshDown)[1])

        mRNAselected_values_newvector <- mRNAselected_values


        if (!is.na(mRNAselected_values_ordered_top)) {
            # How many samples do we have
            numberOfSamples <- length(mRNAselected_values_ordered)

            # High group (above ThreshTop)
            lastelementTOP <-
                max(which(
                    mRNAselected_values_ordered > mRNAselected_values_ordered_top
                ))

            # Low group (below ThreshDown)
            firstelementDOWN <-
                min(
                    which(
                        mRNAselected_values_ordered <= mRNAselected_values_ordered_down
                    )
                )

            samples_top_mRNA_selected <-
                names(mRNAselected_values_ordered[1:lastelementTOP])
            samples_down_mRNA_selected <-
                names(mRNAselected_values_ordered[firstelementDOWN:numberOfSamples])

            # Which samples are in the intermediate group (above ThreshLow and below ThreshTop)
            samples_UNCHANGED_mRNA_selected <-
                names(mRNAselected_values_newvector[which((mRNAselected_values_newvector) > mRNAselected_values_ordered_down &
                                                              mRNAselected_values_newvector < mRNAselected_values_ordered_top
                )])

            cfu_onlyTOP <-
                cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_top_mRNA_selected, ]
            cfu_onlyDOWN <-
                cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_down_mRNA_selected, ]
            cfu_onlyUNCHANGED <-
                cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_UNCHANGED_mRNA_selected, ]

            cfu_ordered <- NULL
            cfu_ordered <- rbind(cfu_onlyTOP, cfu_onlyDOWN)
            cfu <- cfu_ordered

            ttime <- as.numeric(cfu[, "days_to_death"])

            sum(status <- ttime > 0) # morti
            deads_complete <- sum(status <- ttime > 0)

            ttime_only_top <- cfu_onlyTOP[, "days_to_death"]
            deads_top <- sum(ttime_only_top > 0)

            if (dim(cfu_onlyDOWN)[1] >= 1) {
                ttime_only_down <- cfu_onlyDOWN[, "days_to_death"]
                deads_down <- sum(ttime_only_down > 0)
            } else {
                deads_down <- 0
            }

            tabSurv_Matrix[i, "Cancer Deaths"] <- deads_complete
            tabSurv_Matrix[i, "Cancer Deaths with Top"] <- deads_top
            tabSurv_Matrix[i, "Cancer Deaths with Down"] <- deads_down
            tabSurv_Matrix[i, "Mean Normal"] <-
                mean(as.numeric(mRNAselected_values_normal))
            dataCancer_onlyTop_sample <-
                dataCancer[, samples_top_mRNA_selected, drop = FALSE]
            dataCancer_onlyTop_sample_mRNASelected <-
                dataCancer_onlyTop_sample[rownames(dataCancer_onlyTop_sample) == mRNAselected, ]
            dataCancer_onlyDown_sample <-
                dataCancer[, samples_down_mRNA_selected, drop = FALSE]
            dataCancer_onlyDown_sample_mRNASelected <-
                dataCancer_onlyDown_sample[rownames(dataCancer_onlyDown_sample) == mRNAselected, ]
            tabSurv_Matrix[i, "Mean Tumor Top"] <-
                mean(as.numeric(dataCancer_onlyTop_sample_mRNASelected))
            tabSurv_Matrix[i, "Mean Tumor Down"] <-
                mean(as.numeric(dataCancer_onlyDown_sample_mRNASelected))

            ttime[!status] <-
                as.numeric(cfu[!status, "days_to_last_follow_up"])
            ttime[which(ttime == -Inf)] <- 0

            ttime <- survival::Surv(ttime, status)
            rownames(ttime) <- rownames(cfu)
            legendHigh <- paste(mRNAselected, "High")
            legendLow  <- paste(mRNAselected, "Low")
            
            # make a data frame for plotting
            surv_df <- data.frame(ttime = ttime, 
                                  status = status,
                                  exprs_grp = c(rep("low", nrow(cfu_onlyTOP)), rep("high", nrow(cfu_onlyDOWN))))
            tabSurv_pvalue <- tryCatch({
                tabSurv <-
                    survival::survdiff(data = surv_df, ttime  ~ exprs_grp)
                tabSurv_chis <- unlist(tabSurv)$chisq
                tabSurv_pvalue <-
                    as.numeric(1 - pchisq(abs(tabSurv$chisq), df = 1))
            }, error = function(e) {
                return(Inf)
            })
            tabSurv_Matrix[i, "pvalue"] <- tabSurv_pvalue

            if (Survresult == TRUE) {
                # plot survival
                fit <- survival::survfit(data = surv_df, ttime  ~ exprs_grp)
                
                # open new plot
                plot.new() 
                dev.control("enable")
                
                # make plot
                titlePlot <-
                    paste("Kaplan-Meier Survival analysis,\n",
                          cancer_name, "\npvalue=",
                          tabSurv_pvalue)
                plot(
                  fit,
                    col = c("#377eb8", "#e41a1c"),
                    main = titlePlot,
                    xlab = "Days",
                    ylab = "Survival",
                    lwd = 2
                )
                legend(
                    "topright", bty = "n",
                    legend = c(legendLow, legendHigh),
                    col = c("#377eb8", "#e41a1c"),
                    text.col = c("#377eb8", "#e41a1c"),
                    pch = 15
                )
                
                # save the recorded plot
	            	out_plot_list[[mRNAselected]] <- recordPlot()
	            	dev.off()
	            	
            }
        } #end if

    } #end for

    tabSurv_Matrix[tabSurv_Matrix == "-Inf"] <- 0

    tabSurvKM <- tabSurv_Matrix

    # Filtering by selected pvalue < 0.01
    tabSurvKM <- tabSurvKM[tabSurvKM$mRNA != 0, ]
    tabSurvKM <- tabSurvKM[tabSurvKM$pvalue < p.cut, ]
    tabSurvKM <- tabSurvKM[!duplicated(tabSurvKM$mRNA), ]
    rownames(tabSurvKM) <- tabSurvKM$mRNA
    tabSurvKM <- tabSurvKM[, -1]
    tabSurvKM <-
        tabSurvKM[order(tabSurvKM$pvalue, decreasing = FALSE), ]

    colnames(tabSurvKM) <-
        gsub("Cancer", "Group2", colnames(tabSurvKM))
    colnames(tabSurvKM) <-
        gsub("Tumor", "Group2", colnames(tabSurvKM))
    colnames(tabSurvKM) <-
        gsub("Normal", "Group1", colnames(tabSurvKM))
    
    # Replace column names so its more informative
    colnames(tabSurvKM) <- gsub("^Group2 Deaths$", "totalDeath.tumourSamples", colnames(tabSurvKM))
    colnames(tabSurvKM) <- gsub("^Group2 Deaths with Top$", "numDeath.tumourSamples.highExprs", colnames(tabSurvKM))
    colnames(tabSurvKM) <- gsub("^Group2 Deaths with Down$", "numDeath.tumourSamples.lowExprs", colnames(tabSurvKM))
    colnames(tabSurvKM) <- gsub("^Mean Group2 Top$", "meanTpm.tumourSamples.highExprs", colnames(tabSurvKM))
    colnames(tabSurvKM) <- gsub("^Mean Group2 Down$", "meanTpm.tumourSamples.lowExprs", colnames(tabSurvKM))
    colnames(tabSurvKM) <- gsub("^Mean Group1$", "meanTpm.normalSamples", colnames(tabSurvKM))
    
    # Add cancer name to the table
    tabSurvKM$gene_id <- row.names(tabSurvKM)
    tabSurvKM$cancer <- cancer_name
    
    # annotate gene name
    tabSurvKM <- merge(gene_annotation[,c("gene_id", "gene_name")], tabSurvKM, by="gene_id")
    tabSurvKM$p.adjust <- p.adjust(tabSurvKM$pvalue, "BH", n=nrow(tabSurvKM)) 
	
    return(list(stats=tabSurvKM, plots=out_plot_list))
}
