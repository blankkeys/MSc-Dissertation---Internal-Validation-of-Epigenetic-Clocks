# Autosomal probe filtering for GSE87571
# This removes probes annotated outside chromosomes 1-22.
# The autosomal object is saved for beta matrix extraction and modelling.

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

# Keep probes annotated to autosomes only.
autosomal_probes <- annotation$probeID[
  annotation$CpG_chrm %in% paste0("chr", 1:22)
]
mSet_autosomal <- mSet[featureNames(mSet) %in% autosomal_probes, ]

saveRDS(mSet_autosomal, output_mset_file)

dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

autosomal_summary <- data.frame(
  probes_before_autosomal_filter = nrow(mSet),
  autosomal_probes = nrow(mSet_autosomal),
  sex_chromosome_or_unmapped_probes = nrow(mSet) - nrow(mSet_autosomal),
  annotation_source = annotation_file
)

write.csv(
  autosomal_summary,
  summary_file,
  row.names = FALSE
)
