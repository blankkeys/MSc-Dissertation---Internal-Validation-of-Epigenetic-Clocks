# SNP-based probe filtering for GSE51032
# This script removes probes where known SNPs may affect the methylation signal.
# SNPs at or near probe target sites can change probe binding or apparent methylation.
# minfi::dropLociWithSnps() uses the array annotation to remove these affected probes.

library(minfi)

# Load detection-filtered methylation object
mSet <- readRDS("data/GSE51032/mset_normalised_filtered.rds")

probes_before <- nrow(mSet)

# Add genomic annotation so minfi can identify SNP-affected probes.
mSet_genomic <- mapToGenome(mSet)

# Remove probes affected by known SNPs using minfi annotation.
mSet_filtered <- dropLociWithSnps(mSet_genomic)

snp_removed_probes <- setdiff(rownames(mSet), rownames(mSet_filtered))

# Save SNP-filtered methylation object
saveRDS(
  mSet_filtered,
  "data/GSE51032/mset_normalised_filtered_annotation.rds"
)

dir.create("results/external_validation/qc", recursive = TRUE, showWarnings = FALSE)

# Save probes removed by SNP filtering
removed_probes <- data.frame(
  probe_id = snp_removed_probes,
  filter_reason = "SNP_affected_probe"
)

write.csv(
  removed_probes,
  "results/external_validation/qc/gse51032_annotation_removed_probes.csv",
  row.names = FALSE
)

# Save summary
annotation_filtering_summary <- data.frame(
  probes_before_snp_filtering = probes_before,
  snp_affected_probes_removed = length(snp_removed_probes),
  probes_after_snp_filtering = nrow(mSet_filtered),
  snp_filter_rule = "dropLociWithSnps(mSet)",
  filters_not_applied = "cross-reactive, sex chromosome, low-variance"
)

write.csv(
  annotation_filtering_summary,
  "results/external_validation/qc/gse51032_annotation_probe_filtering_summary.csv",
  row.names = FALSE
)

dim(mSet_filtered)
