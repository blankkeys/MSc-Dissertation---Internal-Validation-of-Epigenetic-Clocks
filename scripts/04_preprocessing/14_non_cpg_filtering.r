# Non-CpG probe filtering for GSE87571
# This removes probes with IDs that do not start with cg
# The filtered object is saved for autosomal filtering

library(minfi)

input_mset_file <- "data/GSE87571/mset_normalised_filtered_annotation_crossreactive.rds"
output_mset_file <- "data/GSE87571/mset_normalised_filtered_annotation_crossreactive_cpg.rds"
summary_file <- "results/qc/non_cpg_filtering_summary.csv"

mSet <- readRDS(input_mset_file)

# Keep standard CpG probes only
is_cpg_probe <- grepl("^cg", featureNames(mSet))
mSet_cpg <- mSet[is_cpg_probe, ]

saveRDS(mSet_cpg, output_mset_file)

dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

non_cpg_summary <- data.frame(
  probes_before_filter = nrow(mSet),
  retained_cpg_probes = nrow(mSet_cpg),
  non_cpg_probes_removed = sum(!is_cpg_probe),
  non_cpg_probes_remaining = sum(!grepl("^cg", featureNames(mSet_cpg)))
)

write.csv(
  non_cpg_summary,
  summary_file,
  row.names = FALSE
)
