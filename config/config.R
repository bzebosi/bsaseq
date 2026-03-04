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