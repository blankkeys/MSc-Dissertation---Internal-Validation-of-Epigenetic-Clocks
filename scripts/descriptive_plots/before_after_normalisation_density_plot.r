# Before and after normalisation density plot for GSE87571

library(minfi)

rgSet <- readRDS("data/GSE87571/rgset_raw.rds")
mSet_noob <- readRDS("data/GSE87571/mset_normalised.rds")

# Use 50 samples to keep the plot memory-light
samples_to_plot <- sampleNames(rgSet)[seq_len(min(50, ncol(rgSet)))]

mSet_raw <- preprocessRaw(rgSet[, samples_to_plot])
mSet_noob <- mSet_noob[, samples_to_plot]

dir.create("results/descriptive_plots", recursive = TRUE, showWarnings = FALSE)

pdf("results/descriptive_plots/before_after_normalisation_density_plot.pdf")
par(mfrow = c(1, 2))

densityPlot(
  getBeta(mSet_raw),
  sampGroups = NULL,
  main = "Before normalisation"
)

densityPlot(
  getBeta(mSet_noob),
  sampGroups = NULL,
  main = "After Noob normalisation"
)

dev.off()
