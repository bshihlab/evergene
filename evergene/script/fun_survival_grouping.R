#### Label survival groups

survival_grouping <- function(barcode, exprs, survival_status, survival_time, threshBottom, threshTop){
	# Convert percentage to ratio
	threshTop <- threshTop/100
	threshBottom <- threshBottom/100	
	
	surv_df <- data.frame(barcode = barcode, 
						survival_time = survival_time,
						survival_status = survival_status,
						exprs = exprs)
	# Organise which rows will be kept
	thresholdTop_num <- nrow(surv_df) - floor(nrow(surv_df)*threshTop)
	thresholdBottom_num <- ifelse(threshTop == threshBottom, nrow(surv_df)- thresholdTop_num, floor(nrow(surv_df)*threshBottom))
	thresholdNeither_num <- nrow(surv_df) - thresholdTop_num - thresholdBottom_num
	

	# label top and bottom expression groups			
	surv_df <- surv_df[order(surv_df$exprs, surv_df$barcode),]
	surv_df$exprs_grp <- c(rep("Low.exprs", thresholdBottom_num), rep("Mid.exprs", thresholdNeither_num), rep("High.exprs", thresholdTop_num))
	surv_df$exprs_grp <- factor(surv_df$exprs_grp, levels = c("Low.exprs", "Mid.exprs", "High.exprs"))
	return(surv_df)
}