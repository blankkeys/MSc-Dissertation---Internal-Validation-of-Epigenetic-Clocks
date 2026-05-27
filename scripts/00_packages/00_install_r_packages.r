# Install CRAN packages needed for the analysis
# This installs rsample if it is not already available.

if (!requireNamespace("rsample", quietly = TRUE)) {
  install.packages(
    "rsample",
    repos = "https://cloud.r-project.org"
  )
}
