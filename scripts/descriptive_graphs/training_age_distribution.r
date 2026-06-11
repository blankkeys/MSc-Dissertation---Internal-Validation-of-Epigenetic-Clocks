# Age distribution graph for the GSE87571 training dataset
# This uses the final age-available modelling metadata

metadata <- read.csv("data/GSE87571/modelling_metadata_age_model.csv")

dir.create("results/descriptive_graphs", recursive = TRUE, showWarnings = FALSE)

# Use a histogram to show how many samples fall into each age range
pdf("results/descriptive_graphs/training_age_distribution_histogram.pdf")
hist(
  metadata$age,
  breaks = seq(
    floor(min(metadata$age, na.rm = TRUE) / 10) * 10,
    ceiling(max(metadata$age, na.rm = TRUE) / 10) * 10,
    by = 10
  ),
  xlab = "Chronological age",
  ylab = "Number of samples",
  main = "GSE87571 training dataset age spread"
)
dev.off()

age_distribution_summary <- data.frame(
  samples = nrow(metadata),
  minimum_age = min(metadata$age, na.rm = TRUE),
  first_quartile_age = as.numeric(quantile(metadata$age, 0.25, na.rm = TRUE)),
  median_age = median(metadata$age, na.rm = TRUE),
  mean_age = mean(metadata$age, na.rm = TRUE),
  third_quartile_age = as.numeric(quantile(metadata$age, 0.75, na.rm = TRUE)),
  maximum_age = max(metadata$age, na.rm = TRUE)
)

write.csv(
  age_distribution_summary,
  "results/descriptive_graphs/training_age_distribution_summary.csv",
  row.names = FALSE
)
