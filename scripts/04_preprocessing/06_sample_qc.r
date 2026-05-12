# Sample QC for GSE87571

library(minfi)

# Load the raw RGChannelSet created by 05/05_import_with_minfi.r
rgSet <- readRDS("data/GSE87571/rgset_raw.rds")

# Calculate detection p-values for each probe in each sample.
# A high detection p-value means the probe signal is not clearly above background.
# These values are used to identify samples with many poorly detected probes.
detP <- detectionP(rgSet)

# Save the full detection p-value matrix for later use
saveRDS(detP, "data/GSE87571/sample_detection_pvalues.rds")

# Use p > 0.01 as the initial detection failure threshold.
# For each sample, calculate how many probes failed and what fraction failed.
detection_p_threshold <- 0.01

# create a summary table of sample QC metrics based on detection p-values
sample_qc <- data.frame(
  sample_name = colnames(detP),
  failed_probes = colSums(detP > detection_p_threshold),
  total_probes = nrow(detP)
)

# calculate the fraction of failed probes for each sample,
# a common QC metric to identify poor-quality samples
sample_qc$failed_probe_fraction <- sample_qc$failed_probes / sample_qc$total_probes

# Flag samples where more than 1% of probes failed detection.
# This only flags samples for review; removal happens in the next stage.
sample_qc$sample_qc_status <- ifelse(
  sample_qc$failed_probe_fraction > 0.01,
  "fail",
  "pass"
)

# save the sample QC summary table for later review and filtering
write.csv(
  sample_qc,
  "results/qc/sample_qc_summary.csv",
  row.names = FALSE
)

