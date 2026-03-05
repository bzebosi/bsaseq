############################################################
# 00_setup.R
# Setup environment for BSA-Seq analysis
############################################################

# ---------------------------------------------------------
# Install and load required packages
# ---------------------------------------------------------

source("utilities/install_packages.R")


packages <- c(
 "reshape2", "grid","locfit", "readxl", "BiocManager", "reshape2", "dplyr",
 "zoo", "plyr", "ggplot2","GlobalOptions", "openxlsx", "stringr", 
 "data.table", "naturalsort", "rlang", "purrr", "scales"
)

Install_packages(packages)

source("R/01_read_and_compute.R")
source("R/02_plot_function.R")
source("R/03_visualize_peaks.R")

message("BSA-Seq environment ready.")
