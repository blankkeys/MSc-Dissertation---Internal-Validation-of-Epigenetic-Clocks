# Post-filtering data exploration for GSE87571
# This compares MDS plots before and after sex chromosome probe removal.

library(minfi)


with_sex_chromosomes <- readRDS("data/GSE87571/mset_normalised_filtered_annotation_crossreactive.rds")
without_sex_chromosomes <- readRDS("data/GSE87571/mset_normalised_filtered_annotation_crossreactive_autosomal.rds")


beta_with_sex_chromosomes <- getBeta(with_sex_chromosomes)
beta_without_sex_chromosomes <- getBeta(without_sex_chromosomes)

dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

pdf("results/qc/with_sex_chromosomes_density_plot.pdf")
densityPlot(
  beta_with_sex_chromosomes,
  sampGroups = NULL,
  main = "With sex chromosomes: beta value density"
)
dev.off()

pdf("results/qc/with_sex_chromosomes_mds_plot.pdf")
mdsPlot(
  beta_with_sex_chromosomes,
  numPositions = 10000,
  sampGroups = NULL,
  main = "With sex chromosomes: MDS plot"
)
dev.off()

pdf("results/qc/without_sex_chromosomes_mds_plot.pdf")
mdsPlot(
  beta_without_sex_chromosomes,
  numPositions = 10000,
  sampGroups = NULL,
  main = "Without sex chromosomes: MDS plot"
)
dev.off()

metadata <- read.csv("data/GSE87571/modelling_metadata.csv")
sex_for_samples <- metadata$sex[match(colnames(beta_with_sex_chromosomes), metadata$sample_id)]
samples_with_sex <- !is.na(sex_for_samples)

write.csv(
  as.data.frame(table(metadata$sex, useNA = "ifany")),
  "results/qc/sex_metadata_summary.csv",
  row.names = FALSE
)

pdf("results/qc/with_sex_chromosomes_mds_by_sex.pdf")
mdsPlot(
  beta_with_sex_chromosomes[, samples_with_sex],
  numPositions = 10000,
  sampGroups = factor(sex_for_samples[samples_with_sex]),
  main = "With sex chromosomes: MDS plot by sex"
)
dev.off()

pdf("results/qc/without_sex_chromosomes_mds_by_sex.pdf")
mdsPlot(
  beta_without_sex_chromosomes[, samples_with_sex],
  numPositions = 10000,
  sampGroups = factor(sex_for_samples[samples_with_sex]),
  main = "Without sex chromosomes: MDS plot by sex"
)
dev.off()

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
  "results/qc/post_filtering_data_exploration_summary.csv",
  row.names = FALSE
)
