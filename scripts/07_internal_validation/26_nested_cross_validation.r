# Nested cross-validation for the elastic-net age prediction model
# The outer loop estimates held-out performance.
# The inner loop tunes the elastic-net alpha value.

library(glmnet)
library(rsample)

dir.create("results/internal_validation", recursive = TRUE, showWarnings = FALSE)

set.seed(123)

beta_matrix <- readRDS("data/GSE87571/beta_matrix_age_model.rds")
metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")

# glmnet expects samples as rows and CpG sites as columns.
x <- t(beta_matrix[, match(metadata$sample_id, colnames(beta_matrix))])

# Create nested cross-validation folds.
nested_folds <- nested_cv(
  metadata,
  outside = vfold_cv(v = 10, strata = age),
  inside = vfold_cv(v = 5, strata = age)
)

# Alpha controls the ridge/lasso mixture in glmnet.
# Values from 0.05 to 1.00 test elastic-net mixtures from ridge-like to lasso.
alpha_grid <- seq(0.05, 1, by = 0.05)

outer_performance <- data.frame()
inner_alpha_performance <- data.frame()
outer_residuals <- data.frame()
outer_selected_cpgs <- data.frame()

for (i in seq_len(nrow(nested_folds))) {
  outer_train_metadata <- analysis(nested_folds$splits[[i]])
  outer_test_metadata <- assessment(nested_folds$splits[[i]])
  inner_folds <- nested_folds$inner_resamples[[i]]

  alpha_performance <- data.frame()

  # Test each alpha value within the inner folds.
  for (alpha_value in alpha_grid) {
    for (j in seq_len(nrow(inner_folds))) {
      inner_train_metadata <- analysis(inner_folds$splits[[j]])
      inner_test_metadata <- assessment(inner_folds$splits[[j]])

      x_inner_train <- x[inner_train_metadata$sample_id, ]
      y_inner_train <- inner_train_metadata$age

      x_inner_test <- x[inner_test_metadata$sample_id, ]
      y_inner_test <- inner_test_metadata$age

      inner_model <- cv.glmnet(
        x = x_inner_train,
        y = y_inner_train,
        alpha = alpha_value,
        family = "gaussian"
      )

      inner_predicted_age <- predict(
        inner_model,
        newx = x_inner_test,
        s = "lambda.min"
      )

      inner_predicted_age <- as.numeric(inner_predicted_age)

      fold_alpha_performance <- data.frame(
        outer_fold = nested_folds$id[i],
        inner_fold = inner_folds$id[j],
        alpha = alpha_value,
        lambda_min = inner_model$lambda.min,
        lambda_1se = inner_model$lambda.1se,
        mae = mean(abs(inner_predicted_age - y_inner_test))
      )

      alpha_performance <- rbind(alpha_performance, fold_alpha_performance)
      inner_alpha_performance <- rbind(
        inner_alpha_performance,
        fold_alpha_performance
      )
    }
  }

  # Select the alpha value with the lowest mean inner-fold MAE.
  alpha_summary <- aggregate(
    mae ~ alpha,
    data = alpha_performance,
    FUN = mean
  )

  best_alpha <- alpha_summary$alpha[which.min(alpha_summary$mae)]

  x_outer_train <- x[outer_train_metadata$sample_id, ]
  y_outer_train <- outer_train_metadata$age

  x_outer_test <- x[outer_test_metadata$sample_id, ]
  y_outer_test <- outer_test_metadata$age

  # Train on the full outer training fold using the best alpha from the inner folds.
  outer_model <- cv.glmnet(
    x = x_outer_train,
    y = y_outer_train,
    alpha = best_alpha,
    family = "gaussian"
  )

  # Save the CpGs selected by the outer fold model
  fold_selected_cpgs <- as.matrix(coef(outer_model, s = "lambda.min"))
  fold_selected_cpgs <- data.frame(
    validation_method = "nested_cross_validation",
    resample_id = nested_folds$id[i],
    selected_alpha = best_alpha,
    lambda_min = outer_model$lambda.min,
    lambda_1se = outer_model$lambda.1se,
    cpg = rownames(fold_selected_cpgs),
    coefficient = as.numeric(fold_selected_cpgs[, 1])
  )

  fold_selected_cpgs <- fold_selected_cpgs[
    fold_selected_cpgs$cpg != "(Intercept)" &
      fold_selected_cpgs$coefficient != 0,
  ]

  # Predict the outer held-out fold.
  predicted_age <- predict(
    outer_model,
    newx = x_outer_test,
    s = "lambda.min"
  )

  predicted_age <- as.numeric(predicted_age)

  fold_residuals <- data.frame(
    validation_method = "nested_cross_validation",
    resample_id = nested_folds$id[i],
    sample_id = outer_test_metadata$sample_id,
    geo_accession = outer_test_metadata$geo_accession,
    age = y_outer_test,
    predicted_age = predicted_age,
    residual = predicted_age - y_outer_test,
    absolute_error = abs(predicted_age - y_outer_test)
  )

  fold_performance <- data.frame(
    fold = nested_folds$id[i],
    selected_alpha = best_alpha,
    training_samples = length(y_outer_train),
    test_samples = length(y_outer_test),
    input_cpgs = ncol(x),
    lambda_min = outer_model$lambda.min,
    lambda_1se = outer_model$lambda.1se,
    selected_cpgs = sum(coef(outer_model, s = "lambda.min")[-1, ] != 0),
    mae = mean(abs(predicted_age - y_outer_test)),
    median_absolute_error = median(abs(predicted_age - y_outer_test)),
    rmse = sqrt(mean((predicted_age - y_outer_test)^2)),
    mean_error = mean(predicted_age - y_outer_test),
    correlation = cor(predicted_age, y_outer_test),
    r_squared = cor(predicted_age, y_outer_test)^2
  )

  outer_performance <- rbind(outer_performance, fold_performance)
  outer_residuals <- rbind(outer_residuals, fold_residuals)
  outer_selected_cpgs <- rbind(outer_selected_cpgs, fold_selected_cpgs)
}

nested_cv_summary <- data.frame(
  outer_folds = nrow(outer_performance),
  input_cpgs = ncol(x),
  mean_selected_cpgs = mean(outer_performance$selected_cpgs),
  sd_selected_cpgs = sd(outer_performance$selected_cpgs),
  min_selected_cpgs = min(outer_performance$selected_cpgs),
  max_selected_cpgs = max(outer_performance$selected_cpgs),
  mean_mae = mean(outer_performance$mae),
  sd_mae = sd(outer_performance$mae),
  mean_median_absolute_error = mean(outer_performance$median_absolute_error),
  mean_rmse = mean(outer_performance$rmse),
  sd_rmse = sd(outer_performance$rmse),
  mean_error = mean(outer_performance$mean_error),
  mean_correlation = mean(outer_performance$correlation),
  mean_r_squared = mean(outer_performance$r_squared)
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
