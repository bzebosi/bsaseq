visualize_metrics <- function(
    data, metric = c("af", "ed", "g", "hist", "all"),
    wt = "wildtype", mt = "mutant", prefix = "sample", plots_dir = "plots",
    rollmedian = 501, ylim = NULL, only_mutant = FALSE, device = "png", plot_data = TRUE,
    width = 48, height = 13, hwidth = 30, hheight = 18, dpi = 300, nn_prop = 0.1,
    ems_mode = c("both", "all", "ems"), plot_style = c("grid", "wrap"),
    plot_mode = c("locfit", "rollmedian", "both"), ed_min = 1e-5, g_min = 1e-5, 
    af_min = 0.99, bwidth = 1e6, color_panel = c("blue", "red")
) {
  
  ems_mode <- match.arg(ems_mode)
  plot_style <- match.arg(plot_style)
  plot_mode <- match.arg(plot_mode)
  
  metric <- match.arg(metric, choices = c("af","ed","g","hist","all"), several.ok = TRUE)
  metrics_to_run <- if ("all" %in% metric) c("af", "ed", "g", "hist") else metric
  
  if (length(metrics_to_run) > 1) {
    out <- list()
    for (m in metrics_to_run) {
      out[[m]] <- visualize_metrics(
        data = data, metric = m,
        wt = wt, mt = mt, prefix = prefix, plots_dir = plots_dir,
        rollmedian = rollmedian, ylim = ylim, only_mutant = only_mutant,
        device = device, plot_data = plot_data,
        width = width, height = height, hwidth = hwidth, hheight = hheight,
        dpi = dpi, nn_prop = nn_prop,
        ems_mode = ems_mode, plot_style = plot_style, plot_mode = plot_mode,
        ed_min = ed_min, g_min = g_min, af_min = af_min, bwidth = bwidth,
        color_panel = color_panel
      )
    }
    message("run_metric_only(): multi-metric run completed.")
    return(invisible(out))
  }
  
  # ---- NEW: bmetric must exist (single metric from here on) ----
  bmetric <- metrics_to_run[1]
  
  # pull common tables once
  wt_mt_all  <- data$wt_mt_all
  wt_mt_ems  <- data$wt_mt_ems
  ant_wt_all <- data$ant_wt_all
  ant_wt_ems <- data$ant_wt_ems
  ant_mt_all <- data$ant_mt_all
  ant_mt_ems <- data$ant_mt_ems
  
  # small utilities
  check_cols <- function(df, colname) {
    if (is.null(df)) return(NULL)
    if (!(colname %in% colnames(df))) return(NULL)
    if (nrow(df) == 0) return(NULL)
    df <- df[is.finite(df[[colname]]), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df
  }
  
  filter_min <- function(df, colname, minval = NULL) {
    df <- check_cols(df, colname)
    if (is.null(df)) return(NULL)
    if (!is.null(minval)) {
      df <- df[df[[colname]] >= minval, , drop = FALSE]
      if (nrow(df) == 0) return(NULL)
    }
    df
  }
  
  # shared plotting arguments
  base_args <- list(
    prefix = prefix, plots_dir = plots_dir, ylim = ylim, nn_prop = nn_prop,
    width = width, height = height, hwidth = hwidth, hheight = hheight,
    dpi = dpi, device = device, plot_style = plot_style, color_panel = color_panel
  )
  
  # plot helper
  call_plot <- function(extra_args) {
    args <- modifyList(base_args, extra_args)
    do.call(plot_bsa, args)
  }
  
  # central plot caller (line/point mode)
  call_line <- function(df, column, y_title, plot_title, file_suffix, tr) {
    call_plot(list(
      data = df, column = column, y_title = y_title, plot_title = plot_title,
      file_suffix = file_suffix, is_histogram = FALSE,
      is_locfit = (tr == "locfit"), is_rollmedian = (tr == "rollmedian"),
      rollmedian = if (tr == "rollmedian") rollmedian else 0))
  }
  
  # histogram caller
  call_hist <- function(df, column, y_title, plot_title, file_suffix) {
    call_plot(list(
      data = df, column = column, y_title = y_title, plot_title = plot_title,
      file_suffix = file_suffix, is_histogram = TRUE, threshold = NULL,
      bwidth = bwidth, is_locfit = FALSE, is_rollmedian = FALSE, rollmedian = 0))
  }
  
  af_title  <- function(who, ems = FALSE) sprintf("%s SNP-Index (Unique %sSNPs)", who, if (ems) "ems " else "")
  afd_title <- function(wt, mt) sprintf("AFD: %s - %s", wt, mt)
  who_pair  <- paste0(wt, "_", mt)
  hist_title <- function(who, ems = FALSE) sprintf("%s unique %sSNPs only", who, if (ems) "ems " else "")
  
  
  # parameter builders per metric (returns list of {column, data, plot_title, y_title, suffix_base, snpclass})
  mk <- function(df, column, plot_title, y_title, suffix_base, snpclass, min_filter = NULL) {
    dfx <- if (is.null(min_filter)) check_cols(df, column) else filter_min(df, column, min_filter)
    list(column = column, data = dfx, plot_title = plot_title, y_title = y_title,
         suffix_base = suffix_base, snpclass = snpclass)
  }
  
  make_parameters <- function(bmetric) {
    build_params <- function(ems_flag) {
      wt_mt <- if (ems_flag) wt_mt_ems else wt_mt_all
      wt_df <- if (ems_flag) ant_wt_ems else ant_wt_all
      mt_df <- if (ems_flag) ant_mt_ems else ant_mt_all
      
      snp_shared <- if (ems_flag) "sharedems" else "sharedall"
      snp_uniq   <- if (ems_flag) "uniqems"   else "uniqall"
      
      if (bmetric == "af") {
        if (only_mutant) {
          return(
            list(mk(mt_df, "mt_AF", af_title(mt, ems = ems_flag), "AF",
                    suffix_base = sprintf("AF_%s", mt), snpclass = snp_uniq))
          )
        } else {
          return(
            list(
              mk(wt_mt, "AFD", afd_title(wt, mt), "AFD",
                 suffix_base = sprintf("AFD_%s", who_pair), snpclass = snp_shared),
              mk(wt_df, "wt_AF", af_title(wt, ems = ems_flag), "AF",
                 suffix_base = sprintf("AF_%s", wt), snpclass = snp_uniq),
              mk(mt_df, "mt_AF", af_title(mt, ems = ems_flag), "AF",
                 suffix_base = sprintf("AF_%s", mt), snpclass = snp_uniq))
          )
        }
      }
      
      if (bmetric == "ed") {
        if (only_mutant) {
          message("Skipping ED plots in mutant-only mode (ED requires WT vs MT).")
          return(list())
        } else {
          return(
            list(
              mk(wt_mt, "ED",  sprintf("ED %s & %s",  wt, mt), "ED",
                 suffix_base = sprintf("ED_%s", who_pair), 
                 snpclass = snp_shared, min_filter = ed_min),
              mk(wt_mt, "ED4", sprintf("ED4 %s & %s", wt, mt), "ED4",
                 suffix_base = sprintf("ED4_%s", who_pair), 
                 snpclass = snp_shared, min_filter = ed_min))
          )
        }
      }
      
      if (bmetric == "g") {
        if (only_mutant) {
          message("Skipping G plots in mutant-only mode (G requires WT vs MT).")
          return(list())
        } else {
          return(
            list(
              mk(wt_mt, "G", sprintf("G %s & %s", wt, mt), "G",
                 suffix_base = sprintf("G_%s", who_pair), 
                 snpclass = snp_shared, min_filter = g_min))
          )
        }
      }
      
      if (bmetric == "hist") {
        if (only_mutant) {
          return(
            list(
              mk(mt_df, "mt_AF", hist_title(mt, ems = ems_flag), "SNPs / Mb (×10³)",
                 suffix_base = sprintf("hist_%s", mt), 
                 snpclass = snp_uniq, min_filter = af_min))
          )
        } else {
          return(
            list(
              mk(wt_df, "wt_AF", hist_title(wt, ems = ems_flag), "SNPs / Mb (×10³)",
                 suffix_base = sprintf("hist_%s", wt), snpclass = snp_uniq, min_filter = af_min),
              mk(mt_df, "mt_AF", hist_title(mt, ems = ems_flag), "SNPs / Mb (×10³)",
                 suffix_base = sprintf("hist_%s", mt), snpclass = snp_uniq, min_filter = af_min))
          )
        }
      }
      list()
    }
    
    switch(
      ems_mode, all  = build_params(FALSE), ems  = build_params(TRUE),
      both = c(build_params(FALSE), build_params(TRUE))
    )
  }
  
  parameters <- make_parameters(bmetric)
  
  tracks <- if (plot_mode == "both") c("locfit", "rollmedian") else plot_mode
  
  # plotting loop
  if (plot_data) {
    if (bmetric == "hist") {
      for (pmt in parameters) {
        if (is.null(pmt$data)) next
        message(sprintf("Generating histogram: %s", pmt$plot_title))
        call_hist(
          df = pmt$data, column = pmt$column, y_title = pmt$y_title,plot_title = pmt$plot_title,
          file_suffix = sprintf("%s_%s", pmt$suffix_base, pmt$snpclass)
        )
      }
    } else {
      for (pmt in parameters) {
        if (is.null(pmt$data)) next
        message(sprintf("Generating %s plot: %s", toupper(bmetric), pmt$plot_title))
        for (tr in tracks) {
          if (tr == "rollmedian" && !(rollmedian > 0)) next
          call_line(
            df = pmt$data, column = pmt$column, y_title = pmt$y_title, plot_title = pmt$plot_title,
            file_suffix = sprintf("%s_%s_%s", pmt$suffix_base, tr, pmt$snpclass),
            tr = tr
          )
        }
      }
    }
  }
  invisible(list(parameters = parameters, metric = bmetric, tracks = tracks, ems_mode = ems_mode))
}


run_bsa_one <- function(
    vcf_dir, prefix = "b73", pattern, Genotypes = list(wt = "wildtype", mt = "mutant"),
    min_DP = 5, min_QUAL = 5, only_mutant = FALSE, save_results = FALSE, output_dir = NULL,
    metric = c("af", "ed", "g", "hist", "all"), ems_mode = c("both", "all", "ems"),
    plot_mode = c("locfit", "rollmedian", "both"), plot_style = c("grid", "wrap"),
    plot_data = FALSE, plots_dir = NULL, rollmedian = 501, ylim = NULL, color_panel = c("blue", "red"),
    device = "png", width = 48, height = 13, hwidth = 30, hheight = 18, dpi = 300,
    nn_prop = 0.1, ed_min = 1e-5, g_min = 1e-5, af_min = 0.99, bwidth = 1e6
) {
  
  # ---- ensure dirs exists ----
  ensure_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    invisible(path)
  }
  
  # auto labels from Genotypes
  wt_label <- if (!is.null(Genotypes$wt)) Genotypes$wt else NULL
  mt_label <- Genotypes$mt
  
  # ---- Save or plot ----
  do_save  <- isTRUE(save_results)
  do_plot  <- isTRUE(plot_data)
  
  # ---- match args ----
  ems_mode   <- match.arg(ems_mode)
  plot_mode  <- match.arg(plot_mode)
  plot_style <- match.arg(plot_style)
  metric <- match.arg(metric, choices = c("af", "ed", "g", "hist", "all"), several.ok = TRUE)
  
  # ---- directories (resolve once) ----
  if (is.null(output_dir)) {
    project_dir <- dirname(dirname(vcf_dir))
    output_dir  <- file.path(project_dir, "post_analysis")
  }
  
  # base plots folder
  if (do_save) {
    tab_root <- file.path(output_dir, "tables")
    ensure_dir(tab_root)
    # Create subfolders
    if (only_mutant || is.null(wt_label)) {
      subfolder <- paste(prefix, mt_label, sep = "_")
    } else {
      subfolder <- paste(prefix, wt_label, mt_label, sep = "_")
    }
    tables_dir <- file.path(tab_root, subfolder)
    ensure_dir(tables_dir)
  }
  
  # -- Step 1 : import vcfdata --
  message("=== Step 1: Importing VCF Data === \n ")
  geno_data <- read_vcf(
    vcf_dir = vcf_dir, prefix = prefix, pattern = pattern, Genotypes = Genotypes,
    min_DP = min_DP, min_QUAL = min_QUAL, only_mutant = only_mutant)
  
  # -- Step 1 : compute bsa metrics --
  message("=== Step 2: Analyzing VCF Data === \n")
  res <- compute_bsa(
    geno_data = geno_data, prefix = prefix, save_results = do_save,
    output_dir = if (do_save) tables_dir else output_dir, only_mutant = only_mutant
  )
  
  # -- Step 3 : plot BSA ----
  if (do_plot) {
    # base plots folder
    base_plots_dir <- if (is.null(plots_dir)) file.path(output_dir, "plots") else plots_dir
    ensure_dir(base_plots_dir)
    
    # Create subfolders
    if (only_mutant || is.null(wt_label)) {
      subfolder <- paste(prefix, mt_label, sep = "_")
    } else {
      subfolder <- paste(prefix, wt_label, mt_label, sep = "_")
    }
    plots_dir <- file.path(base_plots_dir, subfolder)
    ensure_dir(plots_dir)
    
    message("=== Step 3: Running BSA Plots === \n")
    res$plot_results <- visualize_metrics(
      data = res, metric = metric, ems_mode = ems_mode, plot_mode = plot_mode,
      wt = wt_label, mt = mt_label, plot_style = plot_style, prefix = prefix, 
      plots_dir = plots_dir, rollmedian = rollmedian, ylim = ylim, 
      only_mutant = only_mutant, device = device, plot_data = TRUE, width = width, 
      height = height, hwidth = hwidth, hheight = hheight, dpi = dpi, nn_prop = nn_prop, 
      ed_min = ed_min, g_min = g_min, af_min = af_min,bwidth = bwidth, color_panel = color_panel
    )
  }
  res
}


run_bsa_all <- function(
    vcf_dir, prefix, pattern, wt_list, mt_list, prefix_list,
    min_DP = 5, min_QUAL = 5, only_mutant = FALSE, save_results = FALSE, output_dir = NULL,
    metric = c("af", "ed", "g", "hist", "all"), ems_mode = c("both", "all", "ems"),
    plot_mode = c("locfit", "rollmedian", "both"), plot_style = c("grid", "wrap"),
    plot_data = FALSE, plots_dir = NULL, rollmedian = 501, ylim = NULL, color_panel = c("blue", "red"),
    device = "png", width = 48, height = 13, hwidth = 30, hheight = 18, dpi = 300,
    nn_prop = 0.1, ed_min = 1e-5, g_min = 1e-5, af_min = 0.99, bwidth = 1e6
) {
  
  results <- list()
  seen_combos <- character()
  
  # ensure output folders exist
  if (is.null(output_dir)) {
    project_dir <- dirname(dirname(vcf_dir))
    output_dir  <- file.path(project_dir, "post_analysis")
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (is.null(plots_dir)) {
    plots_dir <- file.path(output_dir, "plots")
  }
  if (!dir.exists(plots_dir)) dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (pfx in prefix_list) {
    message("\n==== Running BSA pipeline for prefix: ", pfx, " ====")
    if (only_mutant) {
      for (mt in mt_list) {
        combo_key <- paste(pfx, mt, sep = "_")
        if (combo_key %in% seen_combos) next
        seen_combos <- c(seen_combos, combo_key)
        
        message("Running mutant‐only: ", combo_key)
        rex <- run_bsa(
          vcf_dir = vcf_dir, prefix = pfx, pattern = pattern, Genotypes = list(wt = NULL, mt = mt),
          min_DP = min_DP, min_QUAL = min_QUAL, only_mutant = TRUE, 
          save_results = save_results, output_dir = output_dir,metric = metric, 
          ems_mode = ems_mode, plot_mode = plot_mode, plot_style = plot_style,
          plot_data = plot_data, plots_dir = plots_dir, rollmedian = rollmedian, 
          ylim = ylim, color_panel = color_panel,device = device, width = width, 
          height = height, hwidth = hwidth, hheight = hheight, dpi = dpi,
          nn_prop = nn_prop, ed_min = ed_min, g_min = g_min, af_min = af_min, bwidth = bwidth
        )
        results[[combo_key ]] <- rex
      }
    } else {
      # WT vs MT mode: nested loops
      for (wt in wt_list) {
        for (mt in mt_list) {
          if (wt == mt) next
          pair_name <- paste(sort(c(wt, mt)), collapse = "_")
          combo_key <- paste(pfx, pair_name, sep = "_")
          
          if (combo_key %in% seen_combos) next
          seen_combos <- c(seen_combos, combo_key)
          
          message("run:", combo_key, "  (wildtype = ", wt, ",  mutant = ", mt, ")\n")
          rex <- run_bsa_one(
            vcf_dir = vcf_dir, prefix = pfx, pattern = pattern, Genotypes = list(wt = wt, mt = mt),
            min_DP = min_DP, min_QUAL = min_QUAL, only_mutant = FALSE, 
            save_results = save_results, output_dir = output_dir,metric = metric, 
            ems_mode = ems_mode, plot_mode = plot_mode, plot_style = plot_style,
            plot_data = plot_data, plots_dir = plots_dir, rollmedian = rollmedian, 
            ylim = ylim, color_panel = color_panel,device = device, width = width, 
            height = height, hwidth = hwidth, hheight = hheight, dpi = dpi,
            nn_prop = nn_prop, ed_min = ed_min, g_min = g_min, af_min = af_min, bwidth = bwidth
          )
          results[[combo_key ]] <- rex
        }
      }
    }
  }
  message("\n==== run_bsa_all finished ====")
  return(invisible(results))
}