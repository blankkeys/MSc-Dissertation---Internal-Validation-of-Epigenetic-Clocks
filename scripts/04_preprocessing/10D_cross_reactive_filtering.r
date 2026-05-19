# Cross-reactive probe filtering for GSE87571.
# This uses the Zhou/InfiniumAnnotation HM450 manifest after detection p-value,
# non-CpG, and SNP-based filtering have already been completed.

library(minfi)

mSet <- readRDS("data/GSE87571/mset_normalised_filtered_annotation.rds")

annotation_file <- "data/annotation/HM450.hg19.manifest.tsv.gz"
annotation <- read.delim(annotation_file, stringsAsFactors = FALSE, check.names = FALSE)

# Inspect the annotation columns in the Slurm log.
print(names(annotation))

probe_id_column <- "Probe_ID"
mask_column <- "MASK_general"

if (!probe_id_column %in% names(annotation)) {
  stop("Probe_ID column was not found in the Zhou annotation file.")
}

if (!mask_column %in% names(annotation)) {
  stop("MASK_general column was not found. Check the printed column names before choosing a mask column.")
}

annotation <- annotation[match(featureNames(mSet), annotation[[probe_id_column]]), ]

# MASK_general marks probes recommended for masking/removal in the Zhou annotation.
mask_values <- annotation[[mask_column]]
mask_values <- mask_values %in% c(TRUE, "TRUE", "True", "true", 1, "1")

keep <- !mask_values

mSet_filtered <- mSet[keep, ]

saveRDS(
  mSet_filtered,
  "data/GSE87571/mset_normalised_filtered_annotation_crossreactive.rds"
)

removed_probes <- data.frame(
  probe_id = featureNames(mSet)[!keep],
  filter_reason = "cross_reactive_or_non_unique_mapping",
  annotation_column = mask_column,
  annotation_source = annotation_file
)

write.csv(
  removed_probes,
  "results/qc/cross_reactive_removed_probes.csv",
  row.names = FALSE
)

cross_reactive_summary <- data.frame(
  probes_before_cross_reactive_filtering = nrow(mSet),
  probes_removed = sum(!keep),
  probes_after_cross_reactive_filtering = nrow(mSet_filtered),
  annotation_source = annotation_file,
  rule_used = "Remove probes where MASK_general is TRUE"
)

write.csv(
  cross_reactive_summary,
  "results/qc/cross_reactive_probe_filtering_summary.csv",
  row.names = FALSE
)

dim(mSet_filtered)
