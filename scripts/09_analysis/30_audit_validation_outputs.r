# Audit current modelling, validation and external-validation outputs

dir.create("results/analysis", recursive = TRUE, showWarnings = FALSE)

yes_no <- function(value) {
  if (is.na(value)) {
    return("not_applicable")
  }

  if (file.exists(value)) {
    "yes"
  } else {
    "no"
  }
}

read_log_file <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(character())
  }

  readLines(path, warn = FALSE)
}

latest_log_file <- function(log_pattern, extension) {
  if (!dir.exists("logs")) {
    return(NA_character_)
  }

  # Some long validations can be run either as one job or as split/combine jobs
  # log_pattern can therefore be one pattern or several acceptable patterns
  log_files <- unlist(
    lapply(
      log_pattern,
      function(pattern) {
        list.files(
          "logs",
          pattern = paste0("^", pattern, "_[0-9]+(_[0-9]+)?\\.", extension, "$"),
          full.names = TRUE
        )
      }
    )
  )

  if (length(log_files) == 0) {
    return(NA_character_)
  }

  log_files[which.max(file.info(log_files)$mtime)]
}

classify_log_status <- function(log_pattern, summary_file) {
  latest_error_log <- latest_log_file(log_pattern, "err")
  latest_output_log <- latest_log_file(log_pattern, "out")
  error_lines <- read_log_file(latest_error_log)
  output_lines <- read_log_file(latest_output_log)
  combined_lines <- c(error_lines, output_lines)
  combined_text <- paste(combined_lines, collapse = "\n")

  log_exists <- ifelse(
    !is.na(latest_error_log) || !is.na(latest_output_log),
    "yes",
    "no"
  )

  timed_out <- ifelse(
    grepl("TIME LIMIT|CANCELLED AT .* DUE TO TIME LIMIT", combined_text),
    "yes",
    ifelse(log_exists == "yes", "no", "unclear")
  )

  failed <- grepl(
    "Execution halted|Fatal error|Error:|OUT_OF_ME|oom_kill|Killed",
    combined_text,
    ignore.case = TRUE
  )

  summary_exists <- !is.na(summary_file) && file.exists(summary_file)

  completed_successfully <- "unclear"
  if (timed_out == "yes" || failed) {
    completed_successfully <- "no"
  } else if (log_exists == "yes" && summary_exists) {
    completed_successfully <- "yes"
  }

  needs_rerun <- "no"
  if (!summary_exists || completed_successfully == "no") {
    needs_rerun <- "yes"
  } else if (completed_successfully == "unclear") {
    needs_rerun <- "unclear"
  }

  notes <- character()
  if (timed_out == "yes") {
    notes <- c(notes, "latest log indicates time limit")
  }
  if (failed && timed_out != "yes") {
    notes <- c(notes, "latest log indicates failure")
  }
  if (summary_exists && completed_successfully == "no") {
    notes <- c(notes, "existing summary may be from an earlier run")
  }
  if (length(notes) == 0 && completed_successfully == "yes") {
    notes <- "latest available outputs look complete"
  } else if (length(notes) == 0) {
    notes <- "status unclear from available logs"
  }

  list(
    log_file_exists = log_exists,
    latest_error_log = ifelse(is.na(latest_error_log), "", latest_error_log),
    latest_output_log = ifelse(is.na(latest_output_log), "", latest_output_log),
    completed_successfully = completed_successfully,
    timed_out = timed_out,
    needs_rerun = needs_rerun,
    notes = paste(notes, collapse = "; ")
  )
}

audit_stage <- function(
  method,
  script,
  job_wrapper,
  summary_file,
  residual_or_prediction_file,
  selected_cpg_file,
  alpha_tuning_file,
  log_pattern,
  extra_note = ""
) {
  log_status <- classify_log_status(log_pattern, summary_file)

  data.frame(
    method = method,
    script_exists = yes_no(script),
    job_wrapper_exists = yes_no(job_wrapper),
    summary_file_exists = yes_no(summary_file),
    residual_or_prediction_file_exists = yes_no(residual_or_prediction_file),
    selected_cpg_file_exists = yes_no(selected_cpg_file),
    alpha_tuning_file_exists = yes_no(alpha_tuning_file),
    log_file_exists = log_status$log_file_exists,
    completed_successfully = log_status$completed_successfully,
    timed_out = log_status$timed_out,
    needs_rerun = log_status$needs_rerun,
    latest_error_log = log_status$latest_error_log,
    latest_output_log = log_status$latest_output_log,
    notes = paste(c(extra_note, log_status$notes), collapse = "; "),
    stringsAsFactors = FALSE
  )
}

manifest <- rbind(
  audit_stage(
    method = "full_data_tuned_alpha_model",
    script = "scripts/06_elastic_net_modelling/20_elastic_net_regression.r",
    job_wrapper = "jobs/run_20_elastic_net_regression.sh",
    summary_file = "results/modelling/elastic_net_final_model_hyperparameters.csv",
    residual_or_prediction_file = NA,
    selected_cpg_file = "results/modelling/elastic_net_selected_cpg_coefficients.csv",
    alpha_tuning_file = "results/modelling/elastic_net_alpha_tuning_summary.csv",
    log_pattern = "elastic_net_regression"
  ),
  audit_stage(
    method = "full_data_alpha_tuning_summary",
    script = "scripts/06_elastic_net_modelling/20_elastic_net_regression.r",
    job_wrapper = "jobs/run_20_elastic_net_regression.sh",
    summary_file = "results/modelling/elastic_net_alpha_tuning_summary.csv",
    residual_or_prediction_file = NA,
    selected_cpg_file = NA,
    alpha_tuning_file = "results/modelling/elastic_net_alpha_tuning_summary.csv",
    log_pattern = "elastic_net_regression"
  ),
  audit_stage(
    method = "apparent_performance_of_full_data_tuned_alpha_model",
    script = "scripts/07_internal_validation/21_apparent_performance.r",
    job_wrapper = "jobs/run_21_apparent_performance.sh",
    summary_file = "results/internal_validation/apparent_performance_summary.csv",
    residual_or_prediction_file = "results/internal_validation/apparent_performance_residuals.csv",
    selected_cpg_file = "results/internal_validation/apparent_performance_selected_cpgs.csv",
    alpha_tuning_file = NA,
    log_pattern = "apparent_performance",
    extra_note = "apparent performance is optimistic and not an independent validation method"
  ),
  audit_stage(
    method = "single_train_test_split",
    script = "scripts/07_internal_validation/22_single_train_test_split.r",
    job_wrapper = "jobs/run_22_single_train_test_split.sh",
    summary_file = "results/internal_validation/single_train_test_split_summary.csv",
    residual_or_prediction_file = "results/internal_validation/single_train_test_split_residuals.csv",
    selected_cpg_file = "results/internal_validation/single_train_test_split_selected_cpgs.csv",
    alpha_tuning_file = "results/internal_validation/single_train_test_split_alpha_tuning.csv",
    log_pattern = "train_test_split"
  ),
  audit_stage(
    method = "repeated_train_test_split",
    script = "scripts/07_internal_validation/23_repeated_train_test_split.r",
    job_wrapper = "jobs/run_23_repeated_train_test_split.sh",
    summary_file = "results/internal_validation/repeated_train_test_split_summary.csv",
    residual_or_prediction_file = "results/internal_validation/repeated_train_test_split_residuals.csv",
    selected_cpg_file = "results/internal_validation/repeated_train_test_split_selected_cpgs.csv",
    alpha_tuning_file = "results/internal_validation/repeated_train_test_split_alpha_tuning.csv",
    log_pattern = "repeated_train_test"
  ),
  audit_stage(
    method = "k_fold_cross_validation",
    script = "scripts/07_internal_validation/24_k_fold_cross_validation.r",
    job_wrapper = "jobs/run_24_k_fold_cross_validation.sh",
    summary_file = "results/internal_validation/k_fold_summary.csv",
    residual_or_prediction_file = "results/internal_validation/k_fold_residuals.csv",
    selected_cpg_file = "results/internal_validation/k_fold_selected_cpgs.csv",
    alpha_tuning_file = "results/internal_validation/k_fold_alpha_tuning.csv",
    log_pattern = "k_fold_validation"
  ),
  audit_stage(
    method = "repeated_k_fold_cross_validation",
    script = "scripts/07_internal_validation/25_repeated_k_fold_cross_validation.r",
    job_wrapper = "jobs/run_25_repeated_k_fold_cross_validation.sh",
    summary_file = "results/internal_validation/repeated_k_fold_summary.csv",
    residual_or_prediction_file = "results/internal_validation/repeated_k_fold_residuals.csv",
    selected_cpg_file = "results/internal_validation/repeated_k_fold_selected_cpgs.csv",
    alpha_tuning_file = "results/internal_validation/repeated_k_fold_alpha_tuning.csv",
    log_pattern = c(
      "repeated_k_fold_validation",
      "repeated_k_fold_split",
      "combine_repeated_k_fold"
    )
  ),
  audit_stage(
    method = "nested_cross_validation",
    script = "scripts/07_internal_validation/26_nested_cross_validation.r",
    job_wrapper = "jobs/run_26_nested_cross_validation.sh",
    summary_file = "results/internal_validation/nested_cross_validation_summary.csv",
    residual_or_prediction_file = "results/internal_validation/nested_cross_validation_residuals.csv",
    selected_cpg_file = "results/internal_validation/nested_cross_validation_selected_cpgs.csv",
    alpha_tuning_file = "results/internal_validation/nested_cross_validation_inner_alpha_summary.csv",
    log_pattern = c(
      "nested_validation",
      "nested_validation_split",
      "combine_nested_validation"
    )
  ),
  audit_stage(
    method = "bootstrap_632_validation",
    script = "scripts/07_internal_validation/27_bootstrap_632_validation.r",
    job_wrapper = "jobs/run_27_bootstrap_632_validation.sh",
    summary_file = "results/internal_validation/bootstrap_632_summary.csv",
    residual_or_prediction_file = "results/internal_validation/bootstrap_oob_residuals.csv",
    selected_cpg_file = "results/internal_validation/bootstrap_selected_cpgs.csv",
    alpha_tuning_file = "results/internal_validation/bootstrap_alpha_tuning.csv",
    log_pattern = c(
      "bootstrap_632",
      "bootstrap_632_split",
      "combine_bootstrap_632"
    )
  ),
  audit_stage(
    method = "external_gse42861_validation",
    script = "scripts/08_external_validation/08G_external_validation/46_external_gse42861_validation.r",
    job_wrapper = "jobs/run_46_external_gse42861_validation.sh",
    summary_file = "results/external_validation/gse42861_external_validation_summary.csv",
    residual_or_prediction_file = "results/external_validation/gse42861_external_validation_predictions.csv",
    selected_cpg_file = NA,
    alpha_tuning_file = NA,
    log_pattern = "external_validation_GSE42861",
    extra_note = "external validation is evaluation only"
  )
)

write.csv(
  manifest,
  "results/analysis/validation_output_manifest.csv",
  row.names = FALSE
)

print(manifest)
