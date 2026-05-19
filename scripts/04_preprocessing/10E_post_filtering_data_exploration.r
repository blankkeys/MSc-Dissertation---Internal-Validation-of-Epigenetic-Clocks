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

# Check whether sex metadata is present in the methylation object or sample sheet.
metadata <- as.data.frame(pData(mSet), check.names = FALSE)
metadata$sample_id <- rownames(metadata)
metadata$geo_accession <- sub("_.*$", "", metadata$sample_id)

sex_columns <- names(metadata)[grepl("sex|gender", names(metadata), ignore.case = TRUE)]

# If sex is not stored inside the MethylationSet object, check the matched sample sheet.
if (length(sex_columns) == 0 && file.exists("data/GSE87571/qc_ready_sample_sheet.csv")) {
  sample_sheet <- read.csv("data/GSE87571/qc_ready_sample_sheet.csv", check.names = FALSE)
  sample_sheet_sex_columns <- names(sample_sheet)[grepl("sex|gender", names(sample_sheet), ignore.case = TRUE)]

  if (length(sample_sheet_sex_columns) > 0) {
    sex_column <- sample_sheet_sex_columns[1]
    sex_metadata <- data.frame(
      geo_accession = sample_sheet[["geo accession"]],
      sex = sample_sheet[[sex_column]]
    )
    metadata <- merge(metadata, sex_metadata, by = "geo_accession", all.x = TRUE)
  }
} else if (length(sex_columns) > 0) {
  sex_column <- sex_columns[1]
  metadata$sex <- metadata[[sex_column]]
}

if ("sex" %in% names(metadata)) {
  # Standardise sex labels so counts and colours are consistent.
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

  sex_for_beta <- metadata$sex[match(colnames(beta_values), metadata$sample_id)]
  samples_with_sex <- !is.na(sex_for_beta) & sex_for_beta != ""

  if (length(unique(sex_for_beta[samples_with_sex])) > 1) {
    # Create an MDS plot coloured by sex to check whether sex explains sample clustering.
    beta_values_with_sex <- beta_values[, samples_with_sex]
    sex_factor <- factor(sex_for_beta[samples_with_sex])

    pdf("results/qc/post_filtering_mds_by_sex.pdf")
    mdsPlot(
      beta_values_with_sex,
      numPositions = 10000,
      sampGroups = sex_factor,
      main = "Post-filtering MDS plot by sex"
    )
    dev.off()
  }
} else {
  # If sex metadata is absent, record that clearly rather than guessing sex.
  writeLines(
    "No sex or gender metadata column was found in the methylation object or qc_ready_sample_sheet.csv.",
    "results/qc/sex_metadata_summary.txt"
  )
}

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
