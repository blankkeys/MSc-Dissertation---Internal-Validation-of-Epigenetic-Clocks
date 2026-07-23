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

background_cpgs <- unique(rownames(readRDS(background_file)))

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
