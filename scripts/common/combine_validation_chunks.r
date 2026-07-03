# Helper functions for combining Slurm chunk outputs

read_validation_chunks <- function(pattern, expected_chunks) {
  files <- list.files(
    "results/internal_validation",
    pattern = pattern,
    full.names = TRUE
  )
  files <- sort(files)

  if (length(files) != expected_chunks) {
    stop(
      "Expected ",
      expected_chunks,
      " files for ",
      pattern,
      " but found ",
      length(files)
    )
  }

  do.call(rbind, lapply(files, read.csv))
}

summarise_validation_performance <- function(performance, count_name) {
  summary <- data.frame(
    input_cpgs = performance$input_cpgs[1],
    mean_selected_cpgs = mean(performance$selected_cpgs),
    sd_selected_cpgs = sd(performance$selected_cpgs),
    min_selected_cpgs = min(performance$selected_cpgs),
    max_selected_cpgs = max(performance$selected_cpgs),
    mean_selected_alpha = mean(performance$selected_alpha),
    sd_selected_alpha = sd(performance$selected_alpha),
    min_selected_alpha = min(performance$selected_alpha),
    max_selected_alpha = max(performance$selected_alpha),
    mean_lambda_min = mean(performance$lambda_min),
    mean_lambda_1se = mean(performance$lambda_1se),
    mean_mae = mean(performance$mae),
    sd_mae = sd(performance$mae),
    mean_median_absolute_error = mean(performance$median_absolute_error),
    mean_rmse = mean(performance$rmse),
    sd_rmse = sd(performance$rmse),
    mean_error = mean(performance$mean_error),
    mean_correlation = mean(performance$correlation),
    mean_r_squared = mean(performance$r_squared)
  )

  summary <- cbind(
    data.frame(count = nrow(performance)),
    summary
  )
  names(summary)[1] <- count_name

  summary
}
