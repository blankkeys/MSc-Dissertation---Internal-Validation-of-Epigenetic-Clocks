# Combine nested cross-validation chunk outputs

source("scripts/common/combine_validation_chunks.r")

dir.create("results/internal_validation", recursive = TRUE, showWarnings = FALSE)

outer_performance <- read_validation_chunks(
  "^nested_cross_validation_outer_fold_summary_chunk_[0-9]+\\.csv$",
  10
)
outer_residuals <- read_validation_chunks(
  "^nested_cross_validation_residuals_chunk_[0-9]+\\.csv$",
  10
)
inner_alpha_performance <- read_validation_chunks(
  "^nested_cross_validation_inner_alpha_summary_chunk_[0-9]+\\.csv$",
  10
)
outer_selected_cpgs <- read_validation_chunks(
  "^nested_cross_validation_selected_cpgs_chunk_[0-9]+\\.csv$",
  10
)

nested_cv_summary <- summarise_validation_performance(
  outer_performance,
  "outer_folds"
)

write.csv(
  outer_performance,
  "results/internal_validation/nested_cross_validation_outer_fold_summary.csv",
  row.names = FALSE
)

write.csv(
  outer_residuals,
  "results/internal_validation/nested_cross_validation_residuals.csv",
  row.names = FALSE
)

write.csv(
  inner_alpha_performance,
  "results/internal_validation/nested_cross_validation_inner_alpha_summary.csv",
  row.names = FALSE
)

write.csv(
  outer_selected_cpgs,
  "results/internal_validation/nested_cross_validation_selected_cpgs.csv",
  row.names = FALSE
)

write.csv(
  nested_cv_summary,
  "results/internal_validation/nested_cross_validation_summary.csv",
  row.names = FALSE
)
