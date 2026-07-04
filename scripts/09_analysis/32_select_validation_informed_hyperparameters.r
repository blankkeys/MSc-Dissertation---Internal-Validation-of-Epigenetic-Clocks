# Select validation-informed alpha and lambda values for final candidate clocks
# Each internal validation method is reduced to one representative parameter set
# Alpha is selected by the most frequent choice, with lower validation error used to break ties
# Lambda is the median lambda_min among models that selected that alpha

dir.create("results/analysis", recursive = TRUE, showWarnings = FALSE)

method_sources <- data.frame(
  validation_method = c(
    "single_train_test_split",
    "repeated_train_test_split",
    "k_fold_cross_validation",
    "repeated_k_fold_cross_validation",
    "nested_cross_validation",
    "bootstrap_oob"
  ),
  performance_file = c(
    "results/internal_validation/single_train_test_split_summary.csv",
    "results/internal_validation/repeated_train_test_split_per_split_summary.csv",
    "results/internal_validation/k_fold_per_fold_summary.csv",
    "results/internal_validation/repeated_k_fold_per_fold_summary.csv",
    "results/internal_validation/nested_cross_validation_outer_fold_summary.csv",
    "results/internal_validation/bootstrap_632_per_resample_summary.csv"
  ),
  summary_file = c(
    "results/internal_validation/single_train_test_split_summary.csv",
    "results/internal_validation/repeated_train_test_split_summary.csv",
    "results/internal_validation/k_fold_summary.csv",
    "results/internal_validation/repeated_k_fold_summary.csv",
    "results/internal_validation/nested_cross_validation_summary.csv",
    "results/internal_validation/bootstrap_632_summary.csv"
  ),
  performance_column = c(
    "mae",
    "mae",
    "mae",
    "mae",
    "mae",
    "oob_mae"
  ),
  internal_estimate_type = c(
    "held_out_test",
    "repeated_held_out_test",
    "k_fold",
    "repeated_k_fold",
    "nested_outer_fold",
    "bootstrap_oob"
  ),
  stringsAsFactors = FALSE
)

read_if_exists <- function(file) {
  if (!file.exists(file)) {
    warning("Missing file: ", file)
    return(NULL)
  }

  read.csv(file, stringsAsFactors = FALSE)
}

summary_value <- function(summary_data, possible_names) {
  available_name <- possible_names[possible_names %in% names(summary_data)][1]

  if (is.na(available_name)) {
    return(NA_real_)
  }

  summary_data[[available_name]][1]
}

select_alpha <- function(performance_data, performance_column) {
  alpha_counts <- aggregate(
    performance_data[[performance_column]],
    by = list(selected_alpha = performance_data$selected_alpha),
    FUN = length
  )
  names(alpha_counts)[2] <- "selection_count"

  alpha_performance <- aggregate(
    performance_data[[performance_column]],
    by = list(selected_alpha = performance_data$selected_alpha),
    FUN = mean
  )
  names(alpha_performance)[2] <- "mean_validation_error"

  alpha_summary <- merge(alpha_counts, alpha_performance, by = "selected_alpha")
  alpha_summary <- alpha_summary[order(
    -alpha_summary$selection_count,
    alpha_summary$mean_validation_error
  ), ]

  alpha_summary$selected_alpha[1]
}

hyperparameter_selection <- data.frame()

for (i in seq_len(nrow(method_sources))) {
  source <- method_sources[i, ]
  performance_data <- read_if_exists(source$performance_file)
  summary_data <- read_if_exists(source$summary_file)

  if (is.null(performance_data) || is.null(summary_data)) {
    next
  }

  if (!"selected_alpha" %in% names(performance_data)) {
    warning("No selected_alpha column in ", source$performance_file)
    next
  }

  internal_estimated_mae <- summary_value(
    summary_data,
    c("mae", "mean_mae", "mean_oob_mae")
  )
  internal_estimated_rmse <- summary_value(
    summary_data,
    c("rmse", "mean_rmse", "mean_oob_rmse")
  )

  if (source$performance_column %in% names(performance_data)) {
    selected_alpha <- select_alpha(performance_data, source$performance_column)
    selected_alpha_mean_error <- mean(
      performance_data[[source$performance_column]][
        performance_data$selected_alpha == selected_alpha
      ]
    )
  } else {
    selected_alpha <- performance_data$selected_alpha[1]
    selected_alpha_mean_error <- internal_estimated_mae
  }

  selected_alpha_rows <- performance_data$selected_alpha == selected_alpha

  hyperparameter_selection <- rbind(
    hyperparameter_selection,
    data.frame(
      validation_method = source$validation_method,
      internal_estimate_type = source$internal_estimate_type,
      selected_alpha = selected_alpha,
      alpha_selection_count = sum(selected_alpha_rows),
      alpha_selection_total = nrow(performance_data),
      selected_alpha_mean_error = selected_alpha_mean_error,
      selected_lambda_min = median(
        performance_data$lambda_min[selected_alpha_rows],
        na.rm = TRUE
      ),
      selected_lambda_1se = median(
        performance_data$lambda_1se[selected_alpha_rows],
        na.rm = TRUE
      ),
      lambda_selection_rule = "median lambda_min among models selecting chosen alpha",
      internal_estimated_mae = internal_estimated_mae,
      internal_estimated_rmse = internal_estimated_rmse,
      used_external_data_for_selection = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

if (nrow(hyperparameter_selection) == 0) {
  stop("No validation-informed hyperparameters could be selected")
}

write.csv(
  hyperparameter_selection,
  "results/analysis/validation_informed_hyperparameter_selection.csv",
  row.names = FALSE
)

print(hyperparameter_selection)
