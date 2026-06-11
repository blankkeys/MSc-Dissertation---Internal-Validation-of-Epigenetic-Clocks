# Post-filtering data exploration for GSE42861
# This compares CpG-only MDS plots before and after autosomal filtering

library(minfi)

# Load the CpG-filtered object that still contains sex chromosome probes
with_sex_chromosomes <- readRDS("data/GSE42861/mset_normalised_filtered_annotation_crossreactive_cpg.rds")

# Load the CpG-filtered object after autosomal filtering
without_sex_chromosomes <- readRDS("data/GSE42861/mset_normalised_filtered_annotation_crossreactive_cpg_autosomal.rds")

# Extract beta values for plotting
beta_with_sex_chromosomes <- getBeta(with_sex_chromosomes)
beta_without_sex_chromosomes <- getBeta(without_sex_chromosomes)

dir.create("results/external_validation/qc", recursive = TRUE, showWarnings = FALSE)

# Plot beta value density before sex chromosome probe removal
pdf("results/external_validation/qc/gse42861_with_sex_chromosomes_density_plot.pdf")
densityPlot(
  beta_with_sex_chromosomes,
  sampGroups = NULL,
  main = "GSE42861 with sex chromosomes: beta value density"
)
dev.off()

# Plot MDS before sex chromosome probe removal
pdf("results/external_validation/qc/gse42861_with_sex_chromosomes_mds_plot.pdf")
mdsPlot(
  beta_with_sex_chromosomes,
  numPositions = 10000,
  sampGroups = NULL,
  main = "GSE42861 with sex chromosomes: MDS plot"
)
dev.off()

# Plot MDS after sex chromosome probe removal
pdf("results/external_validation/qc/gse42861_without_sex_chromosomes_mds_plot.pdf")
mdsPlot(
  beta_without_sex_chromosomes,
  numPositions = 10000,
  sampGroups = NULL,
  main = "GSE42861 without sex chromosomes: MDS plot"
)
dev.off()

# Save sample and probe counts for both plotted datasets
data_exploration_summary <- data.frame(
  dataset = c("with_sex_chromosomes", "without_sex_chromosomes"),
  samples = c(
    ncol(beta_with_sex_chromosomes),
    ncol(beta_without_sex_chromosomes)
  ),
  probes = c(
    nrow(beta_with_sex_chromosomes),
    nrow(beta_without_sex_chromosomes)
  )
)

write.csv(
  data_exploration_summary,
  "results/external_validation/qc/gse42861_post_filtering_data_exploration_summary.csv",
  row.names = FALSE
)
