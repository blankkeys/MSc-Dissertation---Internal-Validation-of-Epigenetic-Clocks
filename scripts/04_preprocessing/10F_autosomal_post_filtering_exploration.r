# Autosomal post-filtering data exploration for GSE87571
# This removes X/Y chromosome probes for exploratory MDS plots only.
# The main filtered methylation object is not overwritten.

library(minfi)

mSet <- readRDS("data/GSE87571/mset_normalised_filtered_annotation_crossreactive.rds")

annotation_file <- "data/annotation/HM450.hg19.manifest.tsv.gz"
annotation <- read.delim(annotation_file, stringsAsFactors = FALSE, check.names = FALSE)

autosomal_probes <- annotation$probeID[annotation$CpG_chrm %in% paste0("chr", 1:22)]
mSet_autosomal <- mSet[featureNames(mSet) %in% autosomal_probes, ]
beta_values <- getBeta(mSet_autosomal)

dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

# Save an autosomal-only MDS plot without using sex as a grouping variable.
pdf("results/qc/post_filtering_autosomal_mds_plot.pdf")
mdsPlot(
  beta_values,
  numPositions = 10000,
  sampGroups = NULL,
  main = "Post-filtering autosomal MDS plot (ungrouped)"
)
dev.off()

# Read sex metadata from the GEO Series Matrix.
series_lines <- readLines(gzfile("data/GSE87571/GSE87571_series_matrix.txt.gz"))

geo_line <- series_lines[grepl("^!Sample_geo_accession", series_lines)][1]
gender_line <- series_lines[grepl("gender:", series_lines, ignore.case = TRUE)][1]

sex_metadata <- data.frame(
  geo_accession = gsub('"', "", regmatches(geo_line, gregexpr('"[^"]+"', geo_line))[[1]]),
  sex = gsub('"', "", regmatches(gender_line, gregexpr('"[^"]+"', gender_line))[[1]])
)

# Standardise sex labels.
sex_metadata$sex <- tolower(trimws(gsub("sex:|gender:", "", sex_metadata$sex, ignore.case = TRUE)))

# Save an autosomal-only MDS plot coloured by sex.
geo_for_beta <- sub("_.*$", "", colnames(beta_values))
sex_for_beta <- sex_metadata$sex[match(geo_for_beta, sex_metadata$geo_accession)]
beta_values_with_sex <- beta_values[, !is.na(sex_for_beta)]
sex_factor <- factor(sex_for_beta[!is.na(sex_for_beta)])

pdf("results/qc/post_filtering_autosomal_mds_by_sex.pdf")
mdsPlot(
  beta_values_with_sex,
  numPositions = 10000,
  sampGroups = sex_factor,
  main = "Post-filtering autosomal MDS plot by sex"
)
dev.off()

autosomal_summary <- data.frame(
  probes_before_autosomal_filter = nrow(mSet),
  autosomal_probes = nrow(mSet_autosomal),
  sex_chromosome_or_unmapped_probes = nrow(mSet) - nrow(mSet_autosomal),
  annotation_source = annotation_file
)

write.csv(
  autosomal_summary,
  "results/qc/autosomal_post_filtering_exploration_summary.csv",
  row.names = FALSE
)

dim(mSet_autosomal)
