# Before and after normalisation density plot for GSE87571

library(minfi)

# Load raw imported IDAT data
rgSet <- readRDS("data/GSE87571/rgset_raw.rds")

# Load the Noob-normalised methylation object
mSet_noob <- readRDS("data/GSE87571/mset_normalised.rds")

# Convert raw intensities into an unnormalised MethylSet for comparison
mSet_raw <- preprocessRaw(rgSet)

# Extract beta values before and after Noob normalisation
beta_raw <- getBeta(mSet_raw)
beta_noob <- getBeta(mSet_noob)

# Save the plot output folder
dir.create("results/descriptive_plots", recursive = TRUE, showWarnings = FALSE)

# Plot beta-value densities before and after normalisation
pdf("results/descriptive_plots/before_after_normalisation_density_plot.pdf")
par(mfrow = c(1, 2))

densityPlot(
  beta_raw,
  sampGroups = NULL,
  main = "Before normalisation"
)

densityPlot(
  beta_noob,
  sampGroups = NULL,
  main = "After Noob normalisation"
)

dev.off()

# Save the matrix dimensions used in the plot
normalisation_density_summary <- data.frame(
  dataset = c("before_normalisation", "after_noob_normalisation"),
  probes = c(nrow(beta_raw), nrow(beta_noob)),
  samples = c(ncol(beta_raw), ncol(beta_noob))
)

write.csv(
  normalisation_density_summary,
  "results/descriptive_plots/before_after_normalisation_density_summary.csv",
  row.names = FALSE
)
