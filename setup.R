# Install required R packages for the project.
required_packages <- c(
  "CrossScreening",
  "survey",
  "boot",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "scales"
)

missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
} else {
  message("All required R packages are already installed.")
}

message("Python enhanced plots require: pandas, numpy, matplotlib.")
message("Install with: pip install pandas numpy matplotlib")
