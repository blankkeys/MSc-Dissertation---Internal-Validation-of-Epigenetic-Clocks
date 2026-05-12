# remove failed samples from sample qc summary for GSE87571

library(minfi)

# Load the raw RGChannelSet and detection p-values created by previous steps
rgSet <- readRDS("data/GSE87571/rgset_raw.rds")
detp <- readRDS("data/GSE87571/sample_detection_pvalues.rds")
sample_qc <- read.csv("results/qc/sample_qc_summary.csv")

# Identify samples that passed QC based on the sample_qc_summary.csv file
samples_to_keep <- sample_qc$sample_name[sample_qc$sample_qc_status == "pass"]

# Subset the RGChannelSet and detection p-values to keep only the samples that passed QC
rgSet_qc_pass <- rgSet[, samples_to_keep]
detp_qc_pass <- detp[, samples_to_keep]

saveRDS(rgSet_qc_pass, "data/GSE87571/rgset_qc_pass.rds")
saveRDS(detp_qc_pass, "data/GSE87571/sample_detection_pvalues_qc_pass.rds")

