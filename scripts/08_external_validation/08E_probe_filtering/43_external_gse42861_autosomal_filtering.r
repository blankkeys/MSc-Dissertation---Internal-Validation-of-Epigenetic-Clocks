# Autosomal probe filtering for GSE42861
# This removes CpG probes annotated outside chromosomes 1-22
# The filtered object is saved for beta matrix extraction and external validation

library(minfi)

input_mset_file <- "data/GSE42861/mset_normalised_filtered_annotation_crossreactive_cpg.rds"
annotation_file <- "data/annotation/HM450.hg19.manifest.tsv.gz"

output_mset_file <- "data/GSE42861/mset_normalised_filtered_annotation_crossreactive_cpg_autosomal.rds"
summary_file <- "results/external_validation/qc/gse42861_autosomal_filtering_summary.csv"

# load annotation file
mSet <- readRDS(input_mset_file)
annotation <- read.delim(
  annotation_file,
  stringsAsFactors = FALSE, # keep text columns as normal text
  check.names = FALSE # keep original column names from annotation file
)

# Match the annotation rows to the current methylation object
probe_annotation <- annotation[
  match(featureNames(mSet), annotation$probeID),
]

# Keep probes annotated to autosomes only
# mark probes found on chromosomes 1-22
# past0 , join text together with no space
# %in% checks if values are found inside another set of values:
# so, is each probes chromosome one of chr to chr22?
is_autosomal <- probe_annotation$CpG_chrm %in% paste0("chr", 1:22)
# treat probes with missing chromosome annotation as not autosomal
is_autosomal[is.na(is_autosomal)] <- FALSE
# keep only autosomal cpg probes
mSet_autosomal <- mSet[is_autosomal, ]

saveRDS(mSet_autosomal, output_mset_file)

dir.create("results/external_validation/qc", recursive = TRUE, showWarnings = FALSE)

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
