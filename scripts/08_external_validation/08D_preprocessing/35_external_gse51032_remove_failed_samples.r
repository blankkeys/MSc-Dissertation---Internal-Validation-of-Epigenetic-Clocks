# remove failed samples from sample qc summary for GSE51032

library(minfi)

# Load the raw RGChannelSet and detection p-values created by previous steps
rgSet <- readRDS("data/GSE51032/rgset_raw.rds")
detp <- readRDS("data/GSE51032/sample_detection_pvalues.rds")
sample_qc <- read.csv("results/external_validation/qc/gse51032_sample_qc_summary.csv")

# Identify samples that passed QC based on the sample_qc_summary.csv file
samples_to_keep <- sample_qc$sample_name[sample_qc$sample_qc_status == "pass"]

if (length(samples_to_keep) == 0) {
  stop("No GSE51032 samples passed sample QC")
}

# Subset the RGChannelSet and detection p-values to keep only the samples that passed QC
rgSet_qc_pass <- rgSet[, samples_to_keep]
detp_qc_pass <- detp[, samples_to_keep]

saveRDS(rgSet_qc_pass, "data/GSE51032/rgset_qc_pass.rds")
saveRDS(detp_qc_pass, "data/GSE51032/sample_detection_pvalues_qc_pass.rds")
