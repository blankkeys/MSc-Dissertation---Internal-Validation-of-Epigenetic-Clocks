# Probe filtering for GSE51032
# Probes are removed if detection p > 0.01 in more than 1% of samples

library(minfi)

mSet <- readRDS("data/GSE51032/mset_normalised.rds")
detP <- readRDS("data/GSE51032/sample_detection_pvalues_qc_pass.rds")

# Ensure probes are in the same order in the methylation object and detection p-value matrix.
detP <- detP[match(featureNames(mSet), rownames(detP)), ]

# Remove probes that fail detection in more than 1% of samples.
detection_p_threshold <- 0.01
max_failed_fraction <- 0.01

# Calculate the fraction of samples that fail detection for each probe
failed_sample_count <- rowSums(detP > detection_p_threshold)
failed_sample_fraction <- failed_sample_count / ncol(detP)

# Determine which probes to keep based on the failed sample fraction
keep <- failed_sample_fraction <= max_failed_fraction

# Filter the MethylationSet object to retain only the probes that pass the detection p-value filtering
mSet_filtered <- mSet[keep, ]

saveRDS(
  mSet_filtered,
  "data/GSE51032/mset_normalised_filtered.rds"
)

dir.create("results/external_validation/qc", recursive = TRUE, showWarnings = FALSE)

probe_filtering_summary <- data.frame(
  total_probes_before_filtering = nrow(mSet),
  retained_probes = sum(keep),
  removed_probes = sum(!keep),
  detection_p_threshold = detection_p_threshold,
  max_allowed_failed_fraction = max_failed_fraction,
  sample_count_used = ncol(detP)
)

write.csv(
  probe_filtering_summary,
  "results/external_validation/qc/gse51032_probe_filtering_summary.csv",
  row.names = FALSE
)
