#!/bin/bash

# Run the GSE42861 external validation pipeline after download step 29
# This uses sbatch --wait so each job starts after the previous job finishes

set -euo pipefail

sbatch --wait jobs/run_31_check_external_gse42861_metadata.sh
sbatch --wait jobs/run_32_check_external_gse42861_idat.sh
sbatch --wait jobs/run_33_make_external_gse42861_sample_sheet.sh
sbatch --wait jobs/run_34_import_external_gse42861_with_minfi.sh
sbatch --wait jobs/run_35_external_gse42861_sample_qc.sh
sbatch --wait jobs/run_36_external_gse42861_remove_failed_samples.sh
sbatch --wait jobs/run_37_external_gse42861_normalisation.sh
sbatch --wait jobs/run_38_external_gse42861_data_exploration.sh
sbatch --wait jobs/run_39_external_gse42861_det_p_filtering.sh
sbatch --wait jobs/run_40_external_gse42861_annotation_based_filtering.sh
sbatch --wait jobs/run_41_external_gse42861_cross_reactive_filtering.sh
sbatch --wait jobs/run_42_external_gse42861_non_cpg_filtering.sh
sbatch --wait jobs/run_43_external_gse42861_autosomal_filtering.sh
sbatch --wait jobs/run_44_external_gse42861_post_filtering_data_exploration.sh
sbatch --wait jobs/run_45_external_gse42861_extract_beta_matrix.sh
sbatch --wait jobs/run_46_external_gse42861_validation.sh
