base <- "https://raw.githubusercontent.com/bzebosi/bsaseq/main/"

# load the installer helper (note the capital I)
source(paste0(base, "utilities/Install_packages.R"))

# install/load packages (also fix the weird tidyr quote if you copy from 00_setup.R)
packages <- c(
  "reshape2","grid","locfit","readxl","BiocManager","dplyr","tidyr",
  "zoo","plyr","ggplot2","GlobalOptions","openxlsx","stringr",
  "data.table","naturalsort","rlang","purrr","scales","tidyverse"
)
Install_packages(packages)

# load pipeline functions
source(paste0(base, "R/01_read_and_compute.R"))
source(paste0(base, "R/02_plot_function.R"))
source(paste0(base, "R/03_visualize_peaks.R"))

# config: only source gitHub config if you haven't defined it locally
if (!exists("cfg") || !is.list(cfg)) {
  source(paste0(base, "config/config.R"))
}
if (!exists("common_args") || !is.list(common_args)) {
  source(paste0(base, "config/config.R"))
}




# run
results <- list()

for (s in names(cfg)) {
  message("\n=== Running ", s, " ===")
  
  args <- c(
    list(
      vcf_dir = cfg[[s]]$vcf_dir,
      wt_list = cfg[[s]]$wt_list,
      mt_list = cfg[[s]]$mt_list,
      prefix_list = cfg[[s]]$prefix_list
    ),
    common_args
  )
  
  results[[s]] <- do.call(run_bsa_all, args)
}