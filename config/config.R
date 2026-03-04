# Common args
common_args <- list(
  pattern = "snps\\.tsv$", min_DP = 10, min_QUAL = 10,
  metric = "all", plot_mode = "locfit", ems_mode = "both",
  only_mutant = FALSE, save_results = FALSE, plot_data = TRUE
)

# define all samples
cfg <- list(
  S004 = list(
    vcf_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S004/data/snps",
    wt_list = c("S004A"), mt_list = c("S004B"), prefix_list = c("b73")
  ),
  S674 = list(
    vcf_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S674/data/snps",
    wt_list = c("S674A"), mt_list = c("S674B"), prefix_list  = c("b73")
  ),
  S676 = list(
    vcf_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S676/data/snps",
    wt_list = c("S676A"), mt_list = c("S676B"), prefix_list  = c("b73")
  ),
  S687 = list(
    vcf_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S687/data/snps",
    wt_list = c("S687A"), mt_list = c("S687B"), prefix_list  = c("b73")
  )
)

# loop over the samples and run
results <- list()
for (s in names(cfg)) {
  vcf_dir <- cfg[[s]]$vcf_dir
  wt_list <- cfg[[s]]$wt_list
  mt_list <- cfg[[s]]$mt_list
  prefix_list <- cfg[[s]]$prefix_list
  
  message("\n=== Running ", s, " ===")
  cat("\nRunning:", s, "\n",
      "vcf_dir:", vcf_dir, "\n",
      "WT:", paste(wt_list, collapse=", "), "\n",
      "MT:", paste(mt_list, collapse=", "), "\n")
  # run
  args <- c(list(
    vcf_dir = vcf_dir, wt_list = wt_list, mt_list = mt_list, 
    prefix_list = prefix_list), common_args)
  
  results[[s]] <- do.call(run_bsa_all, args)
}
