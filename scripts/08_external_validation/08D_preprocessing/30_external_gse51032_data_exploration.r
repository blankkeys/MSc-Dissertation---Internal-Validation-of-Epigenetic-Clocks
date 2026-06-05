# post normalisation data exploration for GSE51032

library (minfi)

#load the methylation set from normalisation
mSet <- readRDS("data/GSE51032/mset_normalised.rds")

#extract beta values
beta_values <- getBeta(mSet)

dir.create("results/external_validation/qc", recursive = TRUE, showWarnings = FALSE)

#density plot of beta values
pdf("results/external_validation/qc/gse51032_post_norm_density_plots.pdf")
densityPlot(beta_values,
sampGroups = NULL,  # no grouping variable, so all samples will be plotted together
main = "GSE51032 post-normalisation beta value density"
)
dev.off()

#PCA plot of samples using the most variable CpGs
probe_variance <- apply(beta_values, 1, var, na.rm = TRUE)
top_probes <- names(sort(probe_variance, decreasing = TRUE))[1:10000]

pca <- prcomp(t(beta_values[top_probes, ]), scale. = TRUE)

#save PCA plot
pdf("results/external_validation/qc/gse51032_post_norm_pca_plot.pdf")
# Plot the first two principal components
plot(pca$x[, 1],  #first pricipal component, the biggest pattern of variation between samples
pca$x[, 2],  # second principal component, the second biggest variation between sampls
xlab = "PC1",
ylab = "PC2",
main = "GSE51032 post-normalisation PCA",
pch = 16
)
dev.off()   # close the PDF device

#summary statistics of beta values
post_norm_summary <- data.frame(
  samples = ncol(beta_values),
  probes = nrow(beta_values),
  beta_min = min(beta_values, na.rm = TRUE),
  beta_median = median(beta_values, na.rm = TRUE),
  beta_mean = mean(beta_values, na.rm = TRUE),
  beta_max = max(beta_values, na.rm = TRUE)
)

#save summary statistics
write.csv(post_norm_summary, "results/external_validation/qc/gse51032_post_normalisation_summary.csv", row.names = FALSE)
