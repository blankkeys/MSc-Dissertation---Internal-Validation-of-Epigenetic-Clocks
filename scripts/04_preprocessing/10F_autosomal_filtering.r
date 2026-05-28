# Autosomal probe filtering for GSE87571
# This removes CpG probes annotated outside chromosomes 1-22
# The filtered object is saved for beta matrix extraction and modelling

library(minfi)

input_mset_file <- "data/GSE87571/mset_normalised_filtered_annotation_crossreactive_cpg.rds"
annotation_file <- "data/annotation/HM450.hg19.manifest.tsv.gz"

output_mset_file <- "data/GSE87571/mset_normalised_filtered_annotation_crossreactive_cpg_autosomal.rds"
summary_file <- "results/qc/autosomal_filtering_summary.csv"

mSet <- readRDS(input_mset_file)
annotation <- read.delim(
  annotation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Match the annotation rows to the current methylation object
probe_annotation <- annotation[
  match(featureNames(mSet), annotation$probeID),
]

# Keep probes annotated to autosomes only
is_autosomal <- probe_annotation$CpG_chrm %in% paste0("chr", 1:22)
is_autosomal[is.na(is_autosomal)] <- FALSE
mSet_autosomal <- mSet[is_autosomal, ]

saveRDS(mSet_autosomal, output_mset_file)

dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

autosomal_summary <- data.frame(
  probes_before_filter = nrow(mSet),
  retained_autosomal_cpg_probes = nrow(mSet_autosomal),
  sex_chromosome_or_unmapped_probes_removed = sum(!is_autosomal),
  non_cpg_probes_remaining = sum(!grepl("^cg", featureNames(mSet_autosomal))),
  annotation_source = annotation_file
)

write.csv(
  autosomal_summary,
  summary_file,
  row.names = FALSE
)
