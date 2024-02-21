#### Generate example data from gtex
### Data source and preporcessing
# Download count data from https://gtexportal.org/
# Take the first 2000 lines using bash
# head -n 2003 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct > gtex_data_cut.txt
# 

# set working directory
setwd("C:/Users/shihb/OneDrive - Lancaster University/work/project/2023/20230816_nwcr/case_usage/user_input_example")

# read in and select GTEX data
gtex_data <- read.delim("gtex_data_cut.txt")
example_data <- gtex_data[,c("Name", 
                             gsub("-","\\.", c("GTEX-11LCK-0426-SM-5A5M8", "GTEX-11LCK-1226-SM-5Q5AM", "GTEX-11LCK-0926-SM-5A5KA", 
                                               "GTEX-1122O-0126-SM-5GICA", "GTEX-1122O-2426-SM-5GIDN", "GTEX-1122O-0826-SM-5GICV", 
                                               "GTEX-117YX-1326-SM-5H125", "GTEX-117YX-2526-SM-5EQ4Q", "GTEX-117YX-1126-SM-5H128", 
                                               "GTEX-11DXX-0626-SM-5Q5AG", "GTEX-11DXX-2726-SM-5PNXO", "GTEX-11DXX-0326-SM-5PNWC", 
                                               "GTEX-11DXZ-0726-SM-5N9C4", "GTEX-11DXZ-2426-SM-5N9DT", "GTEX-11DXZ-0626-SM-5GU77",
                                               "GTEX-11TT1-1626-SM-5EQL7", "GTEX-11TT1-2326-SM-5GU6N", "GTEX-11TT1-1326-SM-5PNYM",
                                               "GTEX-11GSP-0726-SM-5986L", "GTEX-11GSP-2726-SM-5A5LJ", "GTEX-11GSP-1326-SM-5A5KY")))]
sample_id <- paste0("sample_", 1:21)							
colnames(example_data) <- c("gene_id", sample_id)

subj_id <- paste0("subject", rep(1:7, 3))
sample_annotation <- data.frame(sample_id = sample_id, tissue = rep(c("Lung", "Muscle", "Heart"), 7), subject = subj_id[order(subj_id)])

example_gene_annotation <- gtex_data[,1:2]
colnames(example_gene_annotation) <- c("gene_id", "gene_name")

write.csv(example_data, "example_input.csv", row.names=FALSE)
write.csv(example_gene_annotation, "example_gene_annotation.csv", row.names=FALSE)
write.csv(sample_annotation, "example_sample_annotation.csv", row.names=FALSE)