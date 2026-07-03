# Create a simple age-distribution table and plot for the training dataset
# This is used for presentation slides to show age coverage in GSE87571

dir.create("results/descriptive_plots", recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(
  "data/GSE87571/modelling_metadata_age_model.csv",
  stringsAsFactors = FALSE
)

age_counts <- as.data.frame(table(metadata$age))
names(age_counts) <- c("age", "sample_count")
age_counts$age <- as.numeric(as.character(age_counts$age))

age_summary <- data.frame(
  samples = nrow(metadata),
  minimum_age = min(metadata$age),
  first_quartile_age = as.numeric(quantile(metadata$age, probs = 0.25)),
  median_age = median(metadata$age),
  mean_age = mean(metadata$age),
  third_quartile_age = as.numeric(quantile(metadata$age, probs = 0.75)),
  maximum_age = max(metadata$age)
)

pdf("results/descriptive_plots/training_age_distribution.pdf")
barplot(
  height = age_counts$sample_count,
  names.arg = age_counts$age,
  xlab = "Chronological age",
  ylab = "Number of samples",
  main = "Age distribution of GSE87571 training samples",
  col = "#1f4e79",
  border = NA,
  las = 2,
  cex.names = 0.6
)
dev.off()

write.csv(
  age_counts,
  "results/descriptive_plots/training_age_distribution_counts.csv",
  row.names = FALSE
)

write.csv(
  age_summary,
  "results/descriptive_plots/training_age_distribution_summary.csv",
  row.names = FALSE
)

print(age_summary)
