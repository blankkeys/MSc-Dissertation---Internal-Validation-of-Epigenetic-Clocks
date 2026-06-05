ls #!/bin/bash

# Run the GSE51032 external validation pipeline after download step 22
# This uses sbatch --wait so each job starts after the previous job finishes

set -euo pipefail

sbatch --wait jobs/run_23_check_external_gse51032_metadata.sh
sbatch --wait jobs/run_24_check_external_gse51032_idat.sh
sbatch --wait jobs/run_25_make_external_gse51032_sample_sheet.sh
sbatch --wait jobs/run_26_import_external_gse51032_with_minfi.sh
sbatch --wait jobs/run_27_external_gse51032_sample_qc.sh
sbatch --wait jobs/run_28_external_gse51032_remove_failed_samples.sh
sbatch --wait jobs/run_29_external_gse51032_normalisation.sh
sbatch --wait jobs/run_30_external_gse51032_data_exploration.sh
sbatch --wait jobs/run_31_external_gse51032_det_p_filtering.sh
sbatch --wait jobs/run_32_external_gse51032_annotation_based_filtering.sh
sbatch --wait jobs/run_33_external_gse51032_cross_reactive_filtering.sh
sbatch --wait jobs/run_34_external_gse51032_non_cpg_filtering.sh
sbatch --wait jobs/run_35_external_gse51032_autosomal_filtering.sh
sbatch --wait jobs/run_36_external_gse51032_post_filtering_data_exploration.sh
sbatch --wait jobs/run_37_external_gse51032_extract_beta_matrix.sh
sbatch --wait jobs/run_38_external_gse51032_validation.sh
