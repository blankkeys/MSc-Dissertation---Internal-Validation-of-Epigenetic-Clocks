# Summarise CpG selection stability and alpha/lambda stability
# CpG stability asks whether the same methylation sites are repeatedly selected
# Alpha controls the elastic-net balance between ridge-like and lasso-like penalties
# Lambda controls the overall penalty strength selected by cv.glmnet

normalise_selected_cpgs <- function(selected_cpgs) {
  # Remove intercept rows because the intercept is the model baseline, not a CpG
  selected_cpgs <- selected_cpgs[selected_cpgs$cpg != "(Intercept)", ]

  # Align bootstrap naming with the bootstrap out-of-bag performance summaries
  selected_cpgs$validation_method[
    selected_cpgs$validation_method == "bootstrap"
  ] <- "bootstrap_oob"

  selected_cpgs$coefficient_direction <- ifelse(
    selected_cpgs$coefficient > 0,
    "positive",
    "negative"
  )

  selected_cpgs
}

summarise_method_stability <- function(method_cpgs) {
  # A CpG selected in many folds/resamples is more stable than one selected once
  method <- method_cpgs$validation_method[1]
  model_count <- length(unique(method_cpgs$resample_id))

  cpg_counts <- aggregate(
    method_cpgs$resample_id,
    by = list(cpg = method_cpgs$cpg),
    FUN = length
  )
  names(cpg_counts)[2] <- "selection_count"
  cpg_counts$selection_frequency <- cpg_counts$selection_count / model_count

  data.frame(
    validation_method = method,
    model_count = model_count,
    unique_selected_cpgs = nrow(cpg_counts),
    mean_selected_cpgs_per_model = nrow(method_cpgs) / model_count,
    stable_cpgs_50_percent = sum(cpg_counts$selection_frequency >= 0.50),
    stable_cpgs_80_percent = sum(cpg_counts$selection_frequency >= 0.80),
    stable_cpgs_100_percent = sum(cpg_counts$selection_frequency == 1),
    mean_selection_frequency = mean(cpg_counts$selection_frequency),
    median_selection_frequency = median(cpg_counts$selection_frequency),
    stringsAsFactors = FALSE
  )
}

summarise_cpg_frequency <- function(method_cpgs) {
  # Create a CpG-level table showing how often each CpG was selected
  # This also records whether coefficient direction was stable
  method <- method_cpgs$validation_method[1]
  model_count <- length(unique(method_cpgs$resample_id))

  cpg_counts <- aggregate(
    method_cpgs$resample_id,
    by = list(cpg = method_cpgs$cpg),
    FUN = length
  )
  names(cpg_counts)[2] <- "selection_count"

  mean_coefficients <- aggregate(
    method_cpgs$coefficient,
    by = list(cpg = method_cpgs$cpg),
    FUN = mean
  )
  names(mean_coefficients)[2] <- "mean_coefficient"

  sd_coefficients <- aggregate(
    method_cpgs$coefficient,
    by = list(cpg = method_cpgs$cpg),
    FUN = sd
  )
  names(sd_coefficients)[2] <- "sd_coefficient"

  positive_counts <- aggregate(
    method_cpgs$coefficient > 0,
    by = list(cpg = method_cpgs$cpg),
    FUN = sum
  )
  names(positive_counts)[2] <- "positive_count"

  cpg_frequency <- Reduce(
    function(x, y) merge(x, y, by = "cpg", all = TRUE),
    list(cpg_counts, mean_coefficients, sd_coefficients, positive_counts)
  )
  cpg_frequency$validation_method <- method
  cpg_frequency$model_count <- model_count
  cpg_frequency$selection_frequency <- cpg_frequency$selection_count / model_count
  cpg_frequency$positive_frequency <- cpg_frequency$positive_count /
    cpg_frequency$selection_count
  cpg_frequency$coefficient_direction_stable <- cpg_frequency$positive_frequency == 0 |
    cpg_frequency$positive_frequency == 1

  cpg_frequency[order(-cpg_frequency$selection_frequency), ]
}

calculate_jaccard_similarity <- function(selected_cpgs) {
  # Jaccard similarity measures overlap between two CpG sets
  # Jaccard = shared CpGs divided by total unique CpGs across both sets
  method_sets <- split(selected_cpgs$cpg, selected_cpgs$validation_method)
  method_sets <- lapply(method_sets, unique)
  method_names <- names(method_sets)

  jaccard_similarity <- data.frame()

  for (i in seq_along(method_names)) {
    for (j in seq_along(method_names)) {
      method_a <- method_names[i]
      method_b <- method_names[j]
      shared_cpgs <- length(intersect(method_sets[[method_a]], method_sets[[method_b]]))
      total_cpgs <- length(union(method_sets[[method_a]], method_sets[[method_b]]))

      jaccard_similarity <- rbind(
        jaccard_similarity,
        data.frame(
          method_a = method_a,
          method_b = method_b,
          shared_cpgs = shared_cpgs,
          total_unique_cpgs = total_cpgs,
          jaccard_similarity = shared_cpgs / total_cpgs,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  jaccard_similarity
}

summarise_alpha_lambda <- function(validation_method, file) {
  # Summarise whether each method chooses similar alpha/lambda values across models
  performance_data <- read.csv(file, stringsAsFactors = FALSE)

  if (!"selected_alpha" %in% names(performance_data)) {
    stop("No selected_alpha column in ", file)
  }

  data.frame(
    validation_method = validation_method,
    model_count = nrow(performance_data),
    mean_selected_alpha = mean(performance_data$selected_alpha, na.rm = TRUE),
    sd_selected_alpha = sd(performance_data$selected_alpha, na.rm = TRUE),
    min_selected_alpha = min(performance_data$selected_alpha, na.rm = TRUE),
    max_selected_alpha = max(performance_data$selected_alpha, na.rm = TRUE),
    mean_lambda_min = mean(performance_data$lambda_min, na.rm = TRUE),
    sd_lambda_min = sd(performance_data$lambda_min, na.rm = TRUE),
    mean_lambda_1se = mean(performance_data$lambda_1se, na.rm = TRUE),
    sd_lambda_1se = sd(performance_data$lambda_1se, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

dir.create("results/analysis", recursive = TRUE, showWarnings = FALSE)

# Selected CpG files saved by the model and internal validation scripts
# Repeated methods contain one selected-CpG list per split/fold/resample
selected_cpg_files <- c(
  "results/internal_validation/apparent_performance_selected_cpgs.csv",
  "results/internal_validation/single_train_test_split_selected_cpgs.csv",
  "results/internal_validation/repeated_train_test_split_selected_cpgs.csv",
  "results/internal_validation/k_fold_selected_cpgs.csv",
  "results/internal_validation/repeated_k_fold_selected_cpgs.csv",
  "results/internal_validation/nested_cross_validation_selected_cpgs.csv",
  "results/internal_validation/bootstrap_selected_cpgs.csv"
)

# Performance files used to summarise alpha and lambda stability
# These files contain one row per model, fold, split or bootstrap resample
performance_files <- data.frame(
  validation_method = c(
    "single_train_test_split",
    "repeated_train_test_split",
    "k_fold_cross_validation",
    "repeated_k_fold_cross_validation",
    "nested_cross_validation",
    "bootstrap_oob"
  ),
  file = c(
    "results/internal_validation/single_train_test_split_summary.csv",
    "results/internal_validation/repeated_train_test_split_per_split_summary.csv",
    "results/internal_validation/k_fold_per_fold_summary.csv",
    "results/internal_validation/repeated_k_fold_per_fold_summary.csv",
    "results/internal_validation/nested_cross_validation_outer_fold_summary.csv",
    "results/internal_validation/bootstrap_632_per_resample_summary.csv"
  ),
  stringsAsFactors = FALSE
)

selected_cpgs <- do.call(
  rbind,
  lapply(selected_cpg_files, read.csv, stringsAsFactors = FALSE)
)

selected_cpgs <- normalise_selected_cpgs(selected_cpgs)

cpg_stability_summary <- do.call(
  rbind,
  lapply(
    split(selected_cpgs, selected_cpgs$validation_method),
    summarise_method_stability
  )
)

cpg_frequency_summary <- do.call(
  rbind,
  lapply(
    split(selected_cpgs, selected_cpgs$validation_method),
    summarise_cpg_frequency
  )
)

jaccard_similarity <- calculate_jaccard_similarity(selected_cpgs)

alpha_lambda_stability <- do.call(
  rbind,
  lapply(
    seq_len(nrow(performance_files)),
    function(i) {
      summarise_alpha_lambda(
        performance_files$validation_method[i],
        performance_files$file[i]
      )
    }
  )
)

write.csv(
  cpg_stability_summary,
  "results/analysis/cpg_selection_stability_summary.csv",
  row.names = FALSE
)

write.csv(
  cpg_frequency_summary,
  "results/analysis/cpg_selection_frequency.csv",
  row.names = FALSE
)

write.csv(
  jaccard_similarity,
  "results/analysis/cpg_selection_jaccard_similarity.csv",
  row.names = FALSE
)

write.csv(
  alpha_lambda_stability,
  "results/analysis/alpha_lambda_stability_summary.csv",
  row.names = FALSE
)

print(cpg_stability_summary)
