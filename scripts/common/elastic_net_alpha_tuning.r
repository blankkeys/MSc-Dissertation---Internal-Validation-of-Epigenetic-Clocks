# Helper functions for alpha-tuned elastic-net clock models

tune_alpha_model <- function(x_train, y_train, alpha_grid) {
  nfolds <- min(10, length(y_train))
  foldid <- sample(rep(seq_len(nfolds), length.out = length(y_train)))

  alpha_performance <- data.frame()
  alpha_models <- list()

  for (alpha_value in alpha_grid) {
    elastic_net_model <- cv.glmnet(
      x = x_train,
      y = y_train,
      alpha = alpha_value,
      family = "gaussian",
      foldid = foldid
    )

    alpha_performance <- rbind(
      alpha_performance,
      data.frame(
        alpha = alpha_value,
        lambda_min = elastic_net_model$lambda.min,
        lambda_1se = elastic_net_model$lambda.1se,
        cv_error = min(elastic_net_model$cvm)
      )
    )

    alpha_models[[as.character(alpha_value)]] <- elastic_net_model
  }

  best_alpha <- alpha_performance$alpha[
    which.min(alpha_performance$cv_error)
  ]

  list(
    model = alpha_models[[as.character(best_alpha)]],
    selected_alpha = best_alpha,
    alpha_performance = alpha_performance
  )
}

get_selected_cpgs <- function(model, validation_method, resample_id) {
  selected_cpgs <- as.matrix(coef(model, s = "lambda.min"))
  selected_cpgs <- data.frame(
    validation_method = validation_method,
    resample_id = resample_id,
    cpg = rownames(selected_cpgs),
    coefficient = as.numeric(selected_cpgs[, 1])
  )

  selected_cpgs[
    selected_cpgs$cpg != "(Intercept)" & selected_cpgs$coefficient != 0,
  ]
}
