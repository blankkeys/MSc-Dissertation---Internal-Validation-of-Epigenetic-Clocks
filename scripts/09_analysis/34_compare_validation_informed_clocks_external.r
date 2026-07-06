# Externally validate final clocks in GSE42861 controls
# The benchmark clock uses alpha 0.5 with lambda selected by cv.glmnet
# Validation-parameter clocks use alpha/lambda values carried forward from internal validation
# The external dataset is only used for final testing, not for selecting parameters

dir.create(
  "results/external_validation/validation_informed_clocks",
  recursive = TRUE,
  showWarnings = FALSE
)

# Load external beta values, external metadata and the validation-informed clock coefficients
external_beta_matrix <- readRDS("data/GSE42861/beta_matrix.rds")
training_beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
external_metadata <- read.csv(
  "data/GSE42861/gse42861_qc_ready_sample_sheet_controls.csv",
  stringsAsFactors = FALSE
)
clock_summary <- read.csv(
  "results/modelling/validation_informed_clocks/validation_informed_clock_summary.csv",
  stringsAsFactors = FALSE
)
clock_coefficients <- read.csv(
  "results/modelling/validation_informed_clocks/validation_informed_clock_coefficients.csv",
  stringsAsFactors = FALSE
)

# Choose the metadata sample ID column that matches the beta matrix columns
get_external_sample_id <- function(metadata, beta_matrix) {
  candidate_ids <- list()

  for (column_name in c("Sample_Name", "sample_id", "geo_accession", "Basename")) {
    if (column_name %in% names(metadata)) {
      candidate_ids[[column_name]] <- metadata[[column_name]]
    }
  }

  if ("Basename" %in% names(metadata)) {
    candidate_ids[["Basename_file"]] <- basename(metadata$Basename)
  }

  for (sample_id in candidate_ids) {
    if (any(sample_id %in% colnames(beta_matrix))) {
      return(sample_id)
    }
  }

  stop(
    "No GSE42861 sample ID column matches the beta matrix columns. ",
    "First beta columns: ",
    paste(head(colnames(beta_matrix)), collapse = ", "),
    ". Metadata columns: ",
    paste(names(metadata), collapse = ", ")
  )
}

# Keep only external samples that have both methylation values and age metadata
external_sample_id <- get_external_sample_id(external_metadata, external_beta_matrix)
matched_external_samples <- external_sample_id %in% colnames(external_beta_matrix) &
  !is.na(external_metadata$age)

external_metadata <- external_metadata[matched_external_samples, ]
external_sample_id <- external_sample_id[matched_external_samples]

if (nrow(external_metadata) == 0) {
  stop("No GSE42861 samples with matched beta values and age metadata were found")
}

all_predictions <- data.frame()
external_performance <- data.frame()
all_imputed_predictions <- data.frame()
imputed_external_performance <- data.frame()
compatibility_summary <- data.frame()
missing_cpg_summary <- data.frame()
imputation_summary <- data.frame()

summarise_external_performance <- function(
  method,
  clock_type,
  lambda_source,
  method_clock_summary,
  y_external,
  predicted_age
) {
  AAA <- predicted_age - y_external
  absolute_error <- abs(AAA)
  calibration_model <- lm(predicted_age ~ y_external)
  RAA <- residuals(calibration_model)

  data.frame(
    validation_method = method,
    clock_type = clock_type,
    lambda_source = lambda_source,
    samples = length(y_external),
    selected_alpha = method_clock_summary$selected_alpha[1],
    lambda_min = method_clock_summary$lambda_min[1],
    lambda_1se = method_clock_summary$lambda_1se[1],
    selected_cpgs = method_clock_summary$selected_cpgs[1],
    internal_estimated_mae = method_clock_summary$internal_estimated_mae[1],
    external_mae = mean(absolute_error),
    external_minus_internal_mae = mean(absolute_error) -
      method_clock_summary$internal_estimated_mae[1],
    internal_estimated_rmse = method_clock_summary$internal_estimated_rmse[1],
    external_rmse = sqrt(mean(AAA^2)),
    external_minus_internal_rmse = sqrt(mean(AAA^2)) -
      method_clock_summary$internal_estimated_rmse[1],
    median_absolute_error = median(absolute_error),
    percentile_95_absolute_error = as.numeric(
      quantile(absolute_error, probs = 0.95, na.rm = TRUE)
    ),
    max_absolute_error = max(absolute_error),
    mean_AAA = mean(AAA),
    sd_AAA = sd(AAA),
    mean_RAA = mean(RAA),
    sd_RAA = sd(RAA),
    correlation = cor(predicted_age, y_external),
    r_squared = cor(predicted_age, y_external)^2,
    calibration_intercept = coef(calibration_model)[1],
    calibration_slope = coef(calibration_model)[2],
    stringsAsFactors = FALSE
  )
}

make_external_predictions <- function(
  method,
  clock_type,
  lambda_source,
  external_sample_id,
  external_metadata,
  y_external,
  predicted_age
) {
  AAA <- predicted_age - y_external
  calibration_model <- lm(predicted_age ~ y_external)

  data.frame(
    validation_method = method,
    clock_type = clock_type,
    lambda_source = lambda_source,
    sample_id = external_sample_id,
    geo_accession = external_metadata$geo_accession,
    age = y_external,
    predicted_age = predicted_age,
    residual = AAA,
    AAA = AAA,
    RAA = residuals(calibration_model),
    absolute_error = abs(AAA),
    stringsAsFactors = FALSE
  )
}

summarise_age_bin <- function(age_bin_predictions) {
  data.frame(
    validation_method = age_bin_predictions$validation_method[1],
    age_bin = age_bin_predictions$age_bin[1],
    samples = nrow(age_bin_predictions),
    mean_age = mean(age_bin_predictions$age),
    mae = mean(age_bin_predictions$absolute_error),
    median_absolute_error = median(age_bin_predictions$absolute_error),
    percentile_95_absolute_error = as.numeric(
      quantile(age_bin_predictions$absolute_error, probs = 0.95, na.rm = TRUE)
    ),
    max_absolute_error = max(age_bin_predictions$absolute_error),
    rmse = sqrt(mean(age_bin_predictions$AAA^2)),
    mean_AAA = mean(age_bin_predictions$AAA),
    sd_AAA = sd(age_bin_predictions$AAA),
    mean_RAA = mean(age_bin_predictions$RAA),
    sd_RAA = sd(age_bin_predictions$RAA),
    stringsAsFactors = FALSE
  )
}

# Loop through each fitted clock and apply it to GSE42861
for (method in unique(clock_coefficients$validation_method)) {
  method_coefficients <- clock_coefficients[
    clock_coefficients$validation_method == method,
  ]
  intercept <- method_coefficients$coefficient[
    method_coefficients$cpg == "(Intercept)"
  ]
  selected_coefficients <- method_coefficients[
    method_coefficients$cpg != "(Intercept)",
  ]

  method_clock_summary <- clock_summary[
    clock_summary$validation_method == method,
  ]
  clock_type <- method_clock_summary$clock_type[1]
  lambda_source <- method_clock_summary$lambda_source[1]

# Record whether every selected CpG in the clock exists in the external beta matrix
# If any selected CpG is missing, that specific clock cannot be applied exactly
  missing_cpgs <- setdiff(selected_coefficients$cpg, rownames(external_beta_matrix))
  imputable_cpgs <- intersect(missing_cpgs, rownames(training_beta_matrix))

  compatibility_summary <- rbind(
    compatibility_summary,
    data.frame(
      validation_method = method,
      selected_cpgs = nrow(selected_coefficients),
      available_cpgs = nrow(selected_coefficients) - length(missing_cpgs),
      missing_cpgs = length(missing_cpgs),
      externally_applicable = length(missing_cpgs) == 0,
      mean_imputation_needed = length(missing_cpgs) > 0,
      mean_imputation_applicable = length(missing_cpgs) > 0 &&
        length(missing_cpgs) == length(imputable_cpgs),
      stringsAsFactors = FALSE
    )
  )

  if (length(missing_cpgs) > 0) {
    missing_cpg_summary <- rbind(
      missing_cpg_summary,
      data.frame(
        validation_method = method,
        cpg = missing_cpgs,
        stringsAsFactors = FALSE
      )
    )

    if (length(missing_cpgs) == length(imputable_cpgs)) {
      training_means <- rowMeans(
        training_beta_matrix[missing_cpgs, , drop = FALSE],
        na.rm = TRUE
      )

      imputation_summary <- rbind(
        imputation_summary,
        data.frame(
          validation_method = method,
          cpg = missing_cpgs,
          training_mean_beta = as.numeric(training_means),
          stringsAsFactors = FALSE
        )
      )

      imputed_beta_matrix <- matrix(
        NA_real_,
        nrow = nrow(selected_coefficients),
        ncol = nrow(external_metadata),
        dimnames = list(selected_coefficients$cpg, external_sample_id)
      )

      available_cpgs <- setdiff(selected_coefficients$cpg, missing_cpgs)
      imputed_beta_matrix[available_cpgs, ] <- external_beta_matrix[
        available_cpgs,
        match(external_sample_id, colnames(external_beta_matrix)),
        drop = FALSE
      ]
      imputed_beta_matrix[missing_cpgs, ] <- matrix(
        training_means,
        nrow = length(missing_cpgs),
        ncol = ncol(imputed_beta_matrix)
      )

      x_external_imputed <- t(imputed_beta_matrix)
      y_external <- external_metadata$age
      imputed_predicted_age <- as.numeric(
        intercept + x_external_imputed %*% selected_coefficients$coefficient
      )

      imputed_external_performance <- rbind(
        imputed_external_performance,
        summarise_external_performance(
          method,
          clock_type,
          lambda_source,
          method_clock_summary,
          y_external,
          imputed_predicted_age
        )
      )

      all_imputed_predictions <- rbind(
        all_imputed_predictions,
        make_external_predictions(
          method,
          clock_type,
          lambda_source,
          external_sample_id,
          external_metadata,
          y_external,
          imputed_predicted_age
        )
      )
    }

    next
  }

# Build the external predictor matrix using the selected CpGs in coefficient order
  x_external <- t(external_beta_matrix[
    selected_coefficients$cpg,
    match(external_sample_id, colnames(external_beta_matrix)),
    drop = FALSE
  ])
  y_external <- external_metadata$age

# Predict DNAm age by multiplying external beta values by the trained coefficients
  predicted_age <- as.numeric(intercept + x_external %*% selected_coefficients$coefficient)

# Summarise external prediction performance for this clock
# The internal-external gaps show whether internal validation under- or over-estimated external error
  external_performance <- rbind(
    external_performance,
    summarise_external_performance(
      method,
      clock_type,
      lambda_source,
      method_clock_summary,
      y_external,
      predicted_age
    )
  )

  all_predictions <- rbind(
    all_predictions,
    make_external_predictions(
      method,
      clock_type,
      lambda_source,
      external_sample_id,
      external_metadata,
      y_external,
      predicted_age
    )
  )
}

# Add age bins so performance can be checked across younger and older external samples
if (nrow(all_predictions) > 0) {
  all_predictions$age_bin <- cut(
    all_predictions$age,
    breaks = c(-Inf, 29, 44, 59, 74, Inf),
    labels = c("14-29", "30-44", "45-59", "60-74", "75+"),
    right = TRUE
  )

  external_age_bin_performance <- do.call(
    rbind,
    lapply(
      split(
        all_predictions,
        list(all_predictions$validation_method, all_predictions$age_bin),
        drop = TRUE
      ),
      summarise_age_bin
    )
  )
} else {
  external_age_bin_performance <- data.frame()
}

if (nrow(all_imputed_predictions) > 0) {
  all_imputed_predictions$age_bin <- cut(
    all_imputed_predictions$age,
    breaks = c(-Inf, 29, 44, 59, 74, Inf),
    labels = c("14-29", "30-44", "45-59", "60-74", "75+"),
    right = TRUE
  )

  imputed_external_age_bin_performance <- do.call(
    rbind,
    lapply(
      split(
        all_imputed_predictions,
        list(all_imputed_predictions$validation_method, all_imputed_predictions$age_bin),
        drop = TRUE
      ),
      summarise_age_bin
    )
  )
} else {
  imputed_external_age_bin_performance <- data.frame()
}

# Save CpG compatibility, external predictions and external performance outputs
write.csv(
  compatibility_summary,
  "results/external_validation/validation_informed_clocks/gse42861_clock_compatibility.csv",
  row.names = FALSE
)

write.csv(
  missing_cpg_summary,
  "results/external_validation/validation_informed_clocks/gse42861_missing_clock_cpgs.csv",
  row.names = FALSE
)

write.csv(
  imputation_summary,
  "results/external_validation/validation_informed_clocks/gse42861_training_mean_imputed_cpgs.csv",
  row.names = FALSE
)

write.csv(
  external_performance,
  "results/external_validation/validation_informed_clocks/gse42861_validation_informed_clock_performance.csv",
  row.names = FALSE
)

write.csv(
  all_predictions,
  "results/external_validation/validation_informed_clocks/gse42861_validation_informed_clock_predictions.csv",
  row.names = FALSE
)

write.csv(
  external_age_bin_performance,
  "results/external_validation/validation_informed_clocks/gse42861_validation_informed_clock_age_bin_performance.csv",
  row.names = FALSE
)

write.csv(
  imputed_external_performance,
  "results/external_validation/validation_informed_clocks/gse42861_mean_imputed_clock_performance.csv",
  row.names = FALSE
)

write.csv(
  all_imputed_predictions,
  "results/external_validation/validation_informed_clocks/gse42861_mean_imputed_clock_predictions.csv",
  row.names = FALSE
)

write.csv(
  imputed_external_age_bin_performance,
  "results/external_validation/validation_informed_clocks/gse42861_mean_imputed_clock_age_bin_performance.csv",
  row.names = FALSE
)

print(external_performance)
print(imputed_external_performance)
