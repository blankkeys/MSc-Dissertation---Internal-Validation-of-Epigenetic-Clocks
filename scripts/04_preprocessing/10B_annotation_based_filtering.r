# annotation-based probe filtering for GSE87571
# This script removes non-CpG probes if present and applies SNP-based filtering
# using minfi annotation. Detection p-value filtering was already completed in Stage 10A.

library(minfi)

# Load detection-filtered methylation object
mSet <- readRDS("data/GSE87571/mset_normalised_filtered.rds")

probes_before <- nrow(mSet)

# Remove probes that are not standard CpG probes.
# CpG probe IDs on the 450K array begin with "cg".
non_cpg_probes <- rownames(mSet)[!grepl("^cg", rownames(mSet))]

if (length(non_cpg_probes) > 0) {
  mSet <- mSet[!rownames(mSet) %in% non_cpg_probes, ]
}

probes_after_non_cpg <- nrow(mSet)

# Remove probes with SNPs at the CpG or single-base extension site.
# This uses minfi annotation through dropLociWithSnps().
probes_before_snp <- rownames(mSet)

mSet_snp_filtered <- dropLociWithSnps(
  mSet,
  snps = c("CpG", "SBE"),
  maf = 0
)

snp_removed_probes <- setdiff(probes_before_snp, rownames(mSet_snp_filtered))

mSet <- mSet_snp_filtered

# Save annotation-filtered methylation object
saveRDS(
  mSet,
  "data/GSE87571/mset_normalised_filtered_annotation.rds"
)

# Save probes removed by annotation filtering
removed_probes <- data.frame(
  probe_id = c(non_cpg_probes, snp_removed_probes),
  filter_reason = c(
    rep("non_CpG_probe", length(non_cpg_probes)),
    rep("SNP_at_CpG_or_SBE_site", length(snp_removed_probes))
  )
)

write.csv(
  removed_probes,
  "results/qc/annotation_removed_probes.csv",
  row.names = FALSE
)

# Save summary
annotation_filtering_summary <- data.frame(
  probes_before_annotation_filtering = probes_before,
  non_cpg_probes_removed = length(non_cpg_probes),
  snp_affected_probes_removed = length(snp_removed_probes),
  probes_after_annotation_filtering = nrow(mSet),
  snp_filter_rule = "dropLociWithSnps(snps = c('CpG', 'SBE'), maf = 0)",
  filters_not_applied = "cross-reactive, sex chromosome, low-variance"
)

write.csv(
  annotation_filtering_summary,
  "results/qc/annotation_probe_filtering_summary.csv",
  row.names = FALSE
)

dim(mSet)