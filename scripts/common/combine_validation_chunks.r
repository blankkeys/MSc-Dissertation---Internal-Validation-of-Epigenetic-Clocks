# Helper functions for combining Slurm chunk outputs

# helper function to find files matching a filename pattern and combine them into a single data frame
read_validation_chunks <- function(pattern, expected_chunks) {
  files <- list.files(
    "results/internal_validation",
    pattern = pattern,
    full.names = TRUE
  )
  files <- sort(files)

  # prevent combining incomplete chunks by checking that the number of files matches the expected number of chunks
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

  # read all chunk CSV files and combine them into a single data frame
  chunk_data <- do.call(rbind, lapply(files, read.csv))

  check_alpha_values(chunk_data, "selected_alpha")
  check_alpha_values(chunk_data, "alpha")

  chunk_data
}

# helper function to check that the alpha values in the combined chunks are only those expected 0.25, 0.5, 0.75
# ensure no old alpha values are present
check_alpha_values <- function(chunk_data, column_name) {
  if (!column_name %in% names(chunk_data)) {
    return(invisible(NULL))
  }

  allowed_alpha <- c(0.25, 0.50, 0.75)
  observed_alpha <- sort(unique(chunk_data[[column_name]]))
  unexpected_alpha <- setdiff(observed_alpha, allowed_alpha)

  if (length(unexpected_alpha) > 0) {
    stop(
      "Unexpected alpha values found in combined chunks: ",
      paste(unexpected_alpha, collapse = ", "),
      ". Delete stale chunk files and rerun the split jobs"
    )
  }

  invisible(NULL)
}

# helper function to summarise performance metrics across all chunks
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
  names(summary)[1] <- count_name # turns into a clean summary table by letting the script name the first column properly

  summary
}
