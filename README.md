# Define sample configuration
cfg <- list(
  S004 = list(
    vcf_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S004/data/snps",
    wt_list = c("S004A"),
    mt_list = c("S004B"),
    prefix_list = c("b73")
  )
)

# Define common parameters
common_args <- list(
  pattern = "snps\\.tsv$",
  min_DP = 10,
  min_QUAL = 10,
  metric = "all",
  plot_mode = "locfit",
  ems_mode = "both",
  only_mutant = FALSE,
  save_results = FALSE,
  plot_data = FALSE
)

# Run the pipeline directly from GitHub
source("https://raw.githubusercontent.com/bzebosi/bsaseq/main/scripts/run_bsa_script.R")
