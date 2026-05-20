# Post-filtering data exploration for GSE87571
# This script checks the methylation data after probe filtering.
# It creates density and MDS plots to check the final filtered dataset before modelling.
# If sex metadata is available, it also creates an MDS plot coloured by sex.

library(minfi)

# Load the final filtered MethylationSet object after detection, annotation, and cross-reactive probe filtering.
input_file <- "data/GSE87571/mset_normalised_filtered_annotation_crossreactive.rds"
mSet <- readRDS(input_file)

# Extract beta values for plotting and summary statistics.
beta_values <- getBeta(mSet)

dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

# Save a post-filtering beta value density plot.
pdf("results/qc/post_filtering_density_plot.pdf")
densityPlot(
  beta_values,
  sampGroups = NULL,
  main = "Post-filtering beta value density"
)
dev.off()

# Save an MDS plot, following the Bioconductor methylation workflow style.
pdf("results/qc/post_filtering_mds_plot.pdf")
mdsPlot(
  beta_values,
  numPositions = 10000,
  sampGroups = NULL,
  main = "Post-filtering MDS plot"
)
dev.off()

# Read sex metadata from the GEO Series Matrix.
series_lines <- readLines(gzfile("data/GSE87571/GSE87571_series_matrix.txt.gz"))

geo_line <- series_lines[grepl("^!Sample_geo_accession", series_lines)][1]
gender_line <- series_lines[grepl("gender:", series_lines, ignore.case = TRUE)][1]

geo_accessions <- gsub('"', "", regmatches(geo_line, gregexpr('"[^"]+"', geo_line))[[1]])
sex_values <- gsub('"', "", regmatches(gender_line, gregexpr('"[^"]+"', gender_line))[[1]])

sex_metadata <- data.frame(
  geo_accession = geo_accessions,
  sex = sex_values
)

sample_metadata <- data.frame(
  sample_id = colnames(beta_values),
  geo_accession = sub("_.*$", "", colnames(beta_values))
)

metadata <- merge(sample_metadata, sex_metadata, by = "geo_accession", all.x = TRUE)

# Standardise sex labels.
metadata$sex <- tolower(trimws(gsub("sex:|gender:", "", metadata$sex, ignore.case = TRUE)))
metadata$sex[metadata$sex %in% c("m", "male")] <- "male"
metadata$sex[metadata$sex %in% c("f", "female")] <- "female"

sex_summary <- as.data.frame(table(metadata$sex, useNA = "ifany"))
names(sex_summary) <- c("sex", "sample_count")

write.csv(
  sex_summary,
  "results/qc/sex_metadata_summary.csv",
  row.names = FALSE
)

# Save an MDS plot coloured by sex.
sex_for_beta <- metadata$sex[match(colnames(beta_values), metadata$sample_id)]
beta_values_with_sex <- beta_values[, !is.na(sex_for_beta)]
sex_factor <- factor(sex_for_beta[!is.na(sex_for_beta)])

pdf("results/qc/post_filtering_mds_by_sex.pdf")
mdsPlot(
  beta_values_with_sex,
  numPositions = 10000,
  sampGroups = sex_factor,
  main = "Post-filtering MDS plot by sex"
)
dev.off()

# Save summary statistics for the filtered beta values.
post_filtering_summary <- data.frame(
  input_file = input_file,
  samples = ncol(beta_values),
  probes = nrow(beta_values),
  beta_min = min(beta_values, na.rm = TRUE),
  beta_median = median(beta_values, na.rm = TRUE),
  beta_mean = mean(beta_values, na.rm = TRUE),
  beta_max = max(beta_values, na.rm = TRUE)
)

write.csv(
  post_filtering_summary,
  "results/qc/post_filtering_data_exploration_summary.csv",
  row.names = FALSE
)

dim(mSet)
