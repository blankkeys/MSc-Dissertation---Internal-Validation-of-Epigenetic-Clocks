# Autosomal CpG probe filtering for GSE87571
# This keeps standard CpG probes annotated to chromosomes 1-22
# The filtered object is saved for beta matrix extraction and modelling

library(minfi)

input_mset_file <- "data/GSE87571/mset_normalised_filtered_annotation_crossreactive.rds"
annotation_file <- "data/annotation/HM450.hg19.manifest.tsv.gz"

output_mset_file <- "data/GSE87571/mset_normalised_filtered_annotation_crossreactive_autosomal.rds"
summary_file <- "results/qc/autosomal_filtering_summary.csv"

mSet <- readRDS(input_mset_file)
annotation <- read.delim(
  annotation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Keep standard CpG probes annotated to autosomes only
probe_annotation <- annotation[
  match(featureNames(mSet), annotation$probeID),
]

is_autosomal <- probe_annotation$CpG_chrm %in% paste0("chr", 1:22)
is_cpg_probe <- grepl("^cg", featureNames(mSet))

autosomal_probes <- featureNames(mSet)[
  is_autosomal & is_cpg_probe
]
mSet_autosomal <- mSet[featureNames(mSet) %in% autosomal_probes, ]

saveRDS(mSet_autosomal, output_mset_file)

dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

autosomal_summary <- data.frame(
  probes_before_filter = nrow(mSet),
  retained_autosomal_cpg_probes = nrow(mSet_autosomal),
  sex_chromosome_or_unmapped_probes_removed = sum(!is_autosomal),
  autosomal_non_cpg_probes_removed = sum(is_autosomal & !is_cpg_probe),
  total_removed_probes = nrow(mSet) - nrow(mSet_autosomal),
  non_cpg_probes_remaining = sum(!grepl("^cg", featureNames(mSet_autosomal))),
  annotation_source = annotation_file
)

write.csv(
  autosomal_summary,
  summary_file,
  row.names = FALSE
)
