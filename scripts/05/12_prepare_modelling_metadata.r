# Prepare modelling metadata for GSE87571
# Match age and sex metadata to the final beta matrix samples.

library(GEOquery)

# Load the final beta matrix to get the sample IDs used for modelling.
beta_matrix <- readRDS("data/GSE87571/beta_matrix.rds")

# Extract the GEO sample accession from each beta matrix column name.
modelling_metadata <- data.frame(
  sample_id = colnames(beta_matrix),
  geo_accession = sub("_.*$", "", colnames(beta_matrix))
)

# Read the GEO Series Matrix metadata.
gse <- getGEO(filename = "data/GSE87571/GSE87571_series_matrix.txt.gz")
geo_metadata <- pData(gse)

# Keep the sample accession, age, and sex metadata.
series_metadata <- data.frame(
  geo_accession = rownames(geo_metadata),
  age = geo_metadata$characteristics_ch1,
  sex = geo_metadata$characteristics_ch1.1
)

# Clean age and sex values.
series_metadata$age <- as.numeric(sub("age: ", "", series_metadata$age))
series_metadata$sex <- tolower(sub("gender: ", "", series_metadata$sex))

# Match age and sex to the beta matrix samples using the GEO accession.
metadata_match <- match(modelling_metadata$geo_accession, series_metadata$geo_accession)
modelling_metadata$age <- series_metadata$age[metadata_match]
modelling_metadata$sex <- series_metadata$sex[metadata_match]

# Samples with age available can be used for age prediction modelling.
modelling_metadata$included_in_age_modelling <- !is.na(modelling_metadata$age)

write.csv(
  modelling_metadata,
  "data/GSE87571/modelling_metadata.csv",
  row.names = FALSE
)

dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

# Save a short summary of the modelling metadata.
modelling_metadata_summary <- data.frame(
  beta_matrix_samples = ncol(beta_matrix),
  metadata_rows = nrow(modelling_metadata),
  age_available = sum(!is.na(modelling_metadata$age)),
  age_missing = sum(is.na(modelling_metadata$age)),
  male_samples = sum(modelling_metadata$sex == "male", na.rm = TRUE),
  female_samples = sum(modelling_metadata$sex == "female", na.rm = TRUE),
  samples_included_in_age_modelling = sum(modelling_metadata$included_in_age_modelling),
  minimum_age = min(modelling_metadata$age, na.rm = TRUE),
  maximum_age = max(modelling_metadata$age, na.rm = TRUE)
)

write.csv(
  modelling_metadata_summary,
  "results/qc/modelling_metadata_summary.csv",
  row.names = FALSE
)
