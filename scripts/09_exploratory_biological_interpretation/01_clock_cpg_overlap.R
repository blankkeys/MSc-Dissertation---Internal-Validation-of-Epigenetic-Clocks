# Exploratory CpG overlap analysis with established epigenetic clocks
# This gives biological context only and does not change the main model outputs

known_clock_file <- "data/reference/known_epigenetic_clock_cpgs.csv"
clock_coefficients_file <- "results/modelling/validation_informed_clocks/validation_informed_clock_coefficients.csv"
cpg_frequency_file <- "results/analysis/cpg_selection_frequency.csv"
background_file <- "data/GSE87571/beta_matrix_age_model.rds"
output_dir <- "results/exploratory_cpg_overlap"
stable_80_set <- "stable_cpgs_selected_in_80_percent_models"
stable_100_set <- "stable_cpgs_selected_in_100_percent_models"

dir.create(dirname(known_clock_file), recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(known_clock_file)) {
  write.csv(
    data.frame(clock_name = character(), CpG = character()),
    known_clock_file,
    row.names = FALSE
  )

  stop(
    "Created template at ",
    known_clock_file,
    ". Add published clock CpG lists using columns clock_name,CpG and rerun"
  )
}

known_clock_cpgs <- read.csv(known_clock_file, stringsAsFactors = FALSE)

if (!all(c("clock_name", "CpG") %in% names(known_clock_cpgs))) {
  stop("Known-clock CpG file must contain columns called clock_name and CpG")
}

known_clock_cpgs <- known_clock_cpgs[
  known_clock_cpgs$clock_name != "" & known_clock_cpgs$CpG != "",
]

if (nrow(known_clock_cpgs) == 0) {
  stop("Known-clock CpG file is empty. Add published clock CpG lists and rerun")
}

message("Reading filtered CpG background")
background_cpgs <- unique(rownames(readRDS(background_file)))

message("Reading selected CpGs from validation-informed clocks")
clock_coefficients <- read.csv(clock_coefficients_file, stringsAsFactors = FALSE)
clock_coefficients <- clock_coefficients[clock_coefficients$cpg != "(Intercept)", ]

set_names <- c(
  benchmark_alpha_0.5_cv_lambda = "alpha_0.5_benchmark",
  benchmark_alpha_0.25_cv_lambda = "alpha_0.25_benchmark",
  single_train_test_split = "single_80_20_informed",
  repeated_train_test_split = "repeated_80_20_informed",
  k_fold_cross_validation = "ten_fold_cv_informed",
  repeated_k_fold_cross_validation = "repeated_ten_fold_cv_informed",
  nested_cross_validation = "nested_cv_informed",
  bootstrap_oob = "bootstrap_oob_informed"
)

clock_coefficients$selected_set <- unname(
  set_names[clock_coefficients$validation_method]
)

selected_sets <- split(clock_coefficients$cpg, clock_coefficients$selected_set)
selected_sets <- lapply(selected_sets, unique)

if (file.exists(cpg_frequency_file)) {
  message("Adding stability-defined CpG sets")
  cpg_frequency <- read.csv(cpg_frequency_file, stringsAsFactors = FALSE)

  repeated_methods <- c(
    "repeated_train_test_split",
    "k_fold_cross_validation",
    "repeated_k_fold_cross_validation",
    "nested_cross_validation",
    "bootstrap_oob"
  )

  cpg_frequency <- cpg_frequency[
    cpg_frequency$validation_method %in% repeated_methods,
  ]

  selected_sets[[stable_80_set]] <- unique(
    cpg_frequency$cpg[cpg_frequency$selection_frequency >= 0.80]
  )

  selected_sets[[stable_100_set]] <- unique(
    cpg_frequency$cpg[cpg_frequency$selection_frequency == 1]
  )
}

known_sets <- split(known_clock_cpgs$CpG, known_clock_cpgs$clock_name)
known_sets <- lapply(known_sets, unique)

overlap_summary <- data.frame()

for (selected_set_name in names(selected_sets)) {
  selected_cpgs <- intersect(selected_sets[[selected_set_name]], background_cpgs)

  for (clock_name in names(known_sets)) {
    known_cpgs <- intersect(known_sets[[clock_name]], background_cpgs)
    overlapping_cpgs <- sort(intersect(selected_cpgs, known_cpgs))

    N <- length(background_cpgs)
    K <- length(known_cpgs)
    n <- length(selected_cpgs)
    k <- length(overlapping_cpgs)

    hypergeometric_p_value <- NA
    if (K > 0 && n > 0) {
      hypergeometric_p_value <- phyper(
        k - 1,
        K,
        N - K,
        n,
        lower.tail = FALSE
      )
    }

    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        selected_set = selected_set_name,
        known_clock = clock_name,
        selected_cpgs_in_background = n,
        known_clock_cpgs_in_background = K,
        overlapping_cpgs = k,
        percent_selected_overlapping_known_clock = ifelse(n > 0, 100 * k / n, NA),
        percent_known_clock_recovered = ifelse(K > 0, 100 * k / K, NA),
        jaccard_similarity = ifelse(
          length(union(selected_cpgs, known_cpgs)) > 0,
          k / length(union(selected_cpgs, known_cpgs)),
          NA
        ),
        hypergeometric_p_value = hypergeometric_p_value,
        overlapping_cpg_list = paste(overlapping_cpgs, collapse = ";"),
        stringsAsFactors = FALSE
      )
    )

  }
}

overlap_summary$bh_adjusted_p_value <- p.adjust(
  overlap_summary$hypergeometric_p_value,
  method = "BH"
)

table14_sets <- c(
  "bootstrap_oob_informed",
  "nested_cv_informed",
  "repeated_80_20_informed",
  "repeated_ten_fold_cv_informed",
  "single_80_20_informed",
  "ten_fold_cv_informed",
  stable_80_set,
  stable_100_set
)

table14_overlap_summary <- overlap_summary[
  overlap_summary$selected_set %in% table14_sets,
]

overlap_testing_summary <- data.frame(
  selected_sets_tested = length(selected_sets),
  known_clocks_tested = length(known_sets),
  total_hypergeometric_tests = nrow(overlap_summary),
  table14_selected_sets_shown = length(table14_sets),
  table14_comparisons_shown = nrow(table14_overlap_summary),
  benchmark_comparisons_not_shown = nrow(overlap_summary) -
    nrow(table14_overlap_summary),
  bh_correction_scope = "all hypergeometric tests in clock_cpg_overlap_summary.csv",
  stable_set_rule = paste(
    "union of method-specific CpGs reaching the threshold",
    "within at least one multi-model validation method"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  overlap_summary,
  file.path(output_dir, "clock_cpg_overlap_summary.csv"),
  row.names = FALSE
)

write.csv(
  table14_overlap_summary,
  file.path(output_dir, "clock_cpg_overlap_table14_subset.csv"),
  row.names = FALSE
)

write.csv(
  overlap_testing_summary,
  file.path(output_dir, "clock_cpg_overlap_testing_summary.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "# Exploratory CpG overlap analysis",
    "",
    paste0("Known clock CpGs: ", known_clock_file),
    paste0("Selected clock coefficients: ", clock_coefficients_file),
    paste0("CpG stability file: ", cpg_frequency_file),
    paste0("Filtered background: ", background_file),
    paste0("Background CpGs: ", length(background_cpgs)),
    "",
    paste0("Total selected sets tested: ", length(selected_sets)),
    paste0("Known clocks tested: ", length(known_sets)),
    paste0("Total hypergeometric comparisons: ", nrow(overlap_summary)),
    paste0("Table 14 comparisons shown: ", nrow(table14_overlap_summary)),
    "BH correction was applied once across the complete set of hypergeometric comparisons",
    "",
    "Aggregate stable-set rule:",
    paste(
      "The 80 percent and 100 percent stable sets were created as unions of",
      "method-specific stable CpGs across the five multi-model validation methods"
    ),
    paste0(
      "- 80 percent stable set CpGs: ",
      length(selected_sets[[stable_80_set]])
    ),
    paste0(
      "- 100 percent stable set CpGs: ",
      length(selected_sets[[stable_100_set]])
    ),
    "",
    "Selected sets analysed:",
    paste0("- ", names(selected_sets)),
    "",
    "Known clocks analysed:",
    paste0("- ", names(known_sets)),
    "",
    "Interpretation warning:",
    "Overlap supports biological plausibility but does not prove causal ageing biology",
    "Limited overlap does not invalidate the model because different correlated CpGs can encode similar age-related methylation signal"
  ),
  file.path(output_dir, "exploratory_cpg_overlap_README.md")
)

writeLines(
  c(
    "# Dissertation Text Snippet",
    "",
    "## Methods",
    "Overlap with established epigenetic clocks was assessed as an exploratory biological-context analysis. CpGs selected by the validation-informed clocks and stability-defined CpG sets were compared with published CpG lists from established epigenetic clocks. The filtered autosomal CpG set retained for modelling was used as the background universe because only these CpGs were available for model selection. For the exploratory established-clock overlap analysis, aggregate stable CpG sets were created by taking the union of method-specific stable CpGs across the five multi-model validation procedures. The 80 percent stable set included any CpG selected in at least 80 percent of fitted models within at least one validation method, while the 100 percent stable set included any CpG selected in every fitted model within at least one validation method. For each comparison, overlap count, percentage overlap, Jaccard similarity and hypergeometric enrichment p-value were calculated. The full exploratory analysis evaluated ten selected CpG sets against six established clocks, producing 60 hypergeometric comparisons. Benjamini-Hochberg correction was applied across the complete set of 60 tests.",
    "",
    "## Results Template",
    "The strongest overlap was observed between [SELECTED_SET] and [KNOWN_CLOCK], with [N_OVERLAP] overlapping CpGs. This represented [PERCENT_SELECTED]% of selected CpGs and [PERCENT_CLOCK]% of the published clock CpGs available in the filtered background. All 48 comparisons presented in Table 14 remained statistically significant after Benjamini-Hochberg correction across the complete set of 60 tests.",
    "",
    "## Table Note",
    paste0(
      "Full selected-set ranges summarise results across the six validation-informed final clocks. ",
      "The aggregate 80 percent and 100 percent stable sets were formed by taking the union of CpGs reaching the corresponding within-method selection-frequency threshold in at least one of the five multi-model validation procedures; these sets contained ",
      length(selected_sets[[stable_80_set]]),
      " and ",
      length(selected_sets[[stable_100_set]]),
      " CpGs, respectively. Table 14 presents 48 of the 60 comparisons performed; results for the two benchmark clocks are not shown."
    ),
    "",
    "## Discussion Caution",
    "This analysis is descriptive. CpGs selected by elastic-net clocks are predictive features and should not be interpreted as necessarily causal ageing loci."
  ),
  file.path(output_dir, "dissertation_text_snippet.md")
)

message("CpG overlap analysis complete")
