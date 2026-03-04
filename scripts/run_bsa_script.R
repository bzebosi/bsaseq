source("R/00_setup.R")
source("config/config.R")

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