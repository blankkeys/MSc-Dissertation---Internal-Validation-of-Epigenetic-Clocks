# SNP-based probe filtering for GSE87571
# This script follows the Bioconductor methylation workflow style.
# It uses minfi::dropLociWithSnps() to remove probes affected by known SNPs.
# Detection p-value filtering was already completed in Stage 10A.

library(minfi)

# Load detection-filtered methylation object
mSet <- readRDS("data/GSE87571/mset_normalised_filtered.rds")

probes_before <- nrow(mSet)

# Remove probes affected by known SNPs using minfi annotation.
mSet_filtered <- dropLociWithSnps(mSet)

snp_removed_probes <- setdiff(rownames(mSet), rownames(mSet_filtered))

# Save SNP-filtered methylation object
saveRDS(
  mSet_filtered,
  "data/GSE87571/mset_normalised_filtered_annotation.rds"
)

# Save probes removed by SNP filtering
removed_probes <- data.frame(
  probe_id = snp_removed_probes,
  filter_reason = "SNP_affected_probe"
)

write.csv(
  removed_probes,
  "results/qc/annotation_removed_probes.csv",
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
  "results/qc/annotation_probe_filtering_summary.csv",
  row.names = FALSE
)

dim(mSet_filtered)
