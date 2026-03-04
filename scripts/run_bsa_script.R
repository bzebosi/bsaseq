base <- "https://raw.githubusercontent.com/bzebosi/bsaseq/main/"

# 1) load the installer helper (note the capital I)
source(paste0(base, "utilities/Install_packages.R"))

# 2) install/load packages (also fix the weird tidyr quote if you copy from 00_setup.R)
packages <- c(
  "reshape2","grid","locfit","readxl","BiocManager","dplyr","tidyr",
  "zoo","plyr","ggplot2","GlobalOptions","openxlsx","stringr",
  "data.table","naturalsort","rlang","purrr","scales"
)
Install_packages(packages)

# 3) load pipeline functions
source(paste0(base, "R/01_read_and_compute.R"))
source(paste0(base, "R/02_plot_function.R"))
source(paste0(base, "R/03_visualize_peaks.R"))

# 4) load your config and run the loop (same logic as scripts/run_bsa_script.R)
source(paste0(base, "config/config.R"))

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