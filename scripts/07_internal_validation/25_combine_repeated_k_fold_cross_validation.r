# Combine repeated K-fold cross-validation chunk outputs

source("scripts/common/combine_validation_chunks.r")

dir.create("results/internal_validation", recursive = TRUE, showWarnings = FALSE)

all_performance <- read_validation_chunks(
  "^repeated_k_fold_per_fold_summary_chunk_[0-9]+\\.csv$",
  5
)
all_residuals <- read_validation_chunks(
  "^repeated_k_fold_residuals_chunk_[0-9]+\\.csv$",
  5
)
all_selected_cpgs <- read_validation_chunks(
  "^repeated_k_fold_selected_cpgs_chunk_[0-9]+\\.csv$",
  5
)
all_alpha_tuning <- read_validation_chunks(
  "^repeated_k_fold_alpha_tuning_chunk_[0-9]+\\.csv$",
  5
)

repeated_k_fold_summary <- summarise_validation_performance(
  all_performance,
  "folds"
)

write.csv(
  all_performance,
  "results/internal_validation/repeated_k_fold_per_fold_summary.csv",
  row.names = FALSE
)

write.csv(
  all_residuals,
  "results/internal_validation/repeated_k_fold_residuals.csv",
  row.names = FALSE
)

write.csv(
  all_selected_cpgs,
  "results/internal_validation/repeated_k_fold_selected_cpgs.csv",
  row.names = FALSE
)

write.csv(
  all_alpha_tuning,
  "results/internal_validation/repeated_k_fold_alpha_tuning.csv",
  row.names = FALSE
)

write.csv(
  repeated_k_fold_summary,
  "results/internal_validation/repeated_k_fold_summary.csv",
  row.names = FALSE
)
