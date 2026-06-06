# This script extracts the final beta matrix from the filtered GSE51032 MethylationSet object.
# It checks the matrix is suitable for downstream external validation.

library(minfi)

# Load the filtered MethylationSet object after detection, annotation,
# cross-reactive, autosomal, and standard CpG probe filtering.
filtered_mset_file <- readRDS("data/GSE51032/mset_normalised_filtered_annotation_crossreactive_cpg_autosomal.rds")

# Check the dimensions of the final filtered object.
print(dim(filtered_mset_file))

# Extract the beta matrix from the final filtered object.
beta_matrix <- getBeta(filtered_mset_file)

# Check for missing, duplicated, infinite, or out-of-range values.
missing_value_count <- sum(is.na(beta_matrix))
duplicate_probe_count <- sum(duplicated(rownames(beta_matrix)))
duplicate_sample_count <- sum(duplicated(colnames(beta_matrix)))
infinite_value_count <- sum(is.infinite(beta_matrix))
outside_range_count <- sum(beta_matrix < 0 | beta_matrix > 1, na.rm = TRUE)

min_beta_value <- min(beta_matrix, na.rm = TRUE)
max_beta_value <- max(beta_matrix, na.rm = TRUE)

# Save the beta matrix for external validation.
saveRDS(beta_matrix, "data/GSE51032/beta_matrix.rds")

dir.create("results/external_validation/qc", recursive = TRUE, showWarnings = FALSE)

# Save summary of the final matrix.
final_matrix_summary <- data.frame(
  number_of_probes = nrow(beta_matrix),
  number_of_samples = ncol(beta_matrix),
  missing_value_count = missing_value_count,
  duplicate_probe_id_count = duplicate_probe_count,
  duplicate_sample_id_count = duplicate_sample_count,
  infinite_value_count = infinite_value_count,
  values_outside_0_1_count = outside_range_count,
  minimum_beta_value = min_beta_value,
  maximum_beta_value = max_beta_value
)

write.csv(
  final_matrix_summary,
  "results/external_validation/qc/gse51032_final_matrix_summary.csv",
  row.names = FALSE
)
