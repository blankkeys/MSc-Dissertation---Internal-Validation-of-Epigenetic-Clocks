# Predicted age versus chronological age plots for selected validation methods

dir.create("results/descriptive_plots", recursive = TRUE, showWarnings = FALSE)

plot_predicted_vs_age <- function(predictions, output_file, plot_title) {
  predictions <- na.omit(predictions[c("age", "predicted_age")])
  error <- predictions$predicted_age - predictions$age
  calibration <- lm(predicted_age ~ age, data = predictions)

  plot_range <- range(
    c(predictions$age, predictions$predicted_age),
    na.rm = TRUE
  )

  pdf(output_file, width = 6, height = 6)
  plot(
    predictions$age,
    predictions$predicted_age,
    xlim = plot_range,
    ylim = plot_range,
    xlab = "Chronological age",
    ylab = "Predicted DNAm age",
    main = plot_title,
    pch = 16,
    col = rgb(31, 94, 150, 120, maxColorValue = 255)
  )
  abline(0, 1, col = "red", lwd = 2)
  abline(calibration, col = "black", lwd = 2, lty = 2)

  legend(
    "topleft",
    legend = c(
      sprintf("MAE = %.2f years", mean(abs(error))),
      sprintf("RMSE = %.2f years", sqrt(mean(error^2))),
      sprintf("Pearson r = %.3f", cor(predictions$age, predictions$predicted_age)),
      sprintf("Calibration slope = %.3f", coef(calibration)[2]),
      sprintf("Calibration intercept = %.2f", coef(calibration)[1]),
      paste0("n = ", nrow(predictions))
    ),
    bty = "o",
    bg = "white",
    cex = 0.75
  )

  legend(
    "bottomright",
    legend = c("Perfect agreement", "Fitted calibration"),
    col = c("red", "black"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n",
    cex = 0.75
  )
  dev.off()
}

apparent_predictions <- read.csv(
  "results/internal_validation/apparent_performance_residuals.csv"
)

k_fold_predictions <- read.csv(
  "results/internal_validation/k_fold_residuals.csv"
)

bootstrap_oob_predictions <- read.csv(
  "results/internal_validation/bootstrap_oob_residuals.csv"
)

plot_predicted_vs_age(
  apparent_predictions,
  "results/descriptive_plots/apparent_predicted_vs_chronological_age.pdf",
  "Apparent performance"
)

plot_predicted_vs_age(
  k_fold_predictions,
  "results/descriptive_plots/k_fold_predicted_vs_chronological_age.pdf",
  "10-fold cross-validation"
)

plot_predicted_vs_age(
  bootstrap_oob_predictions,
  "results/descriptive_plots/bootstrap_oob_predicted_vs_chronological_age.pdf",
  "Bootstrap out-of-bag validation"
)
