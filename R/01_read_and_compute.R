#' Read SNP tables and compute bulk allele frequency
#' Reads SNP tables (TSV/CSV) for wild-type and mutant bulks, filters SNPs,
#' computes allele frequency (AF), and returns a list of data.tables.
#' @param vcf_dir Directory containing SNP tables.
#' @param prefix Prefix used to match filenames (e.g., "b73").
#' @param pattern Regex to identify candidate files (e.g., "snps\\.tsv$").
#' @param Genotypes Named list defining labels used in filenames.
#' @param min_DP Minimum read depth (DP) filter.
#' @param min_QUAL Minimum quality (QUAL) filter.
#' @param only_mutant If TRUE, only the mutant bulk is loaded.
#' @return A named list of data.tables (wt and mt, or mt only).
#' @export
read_vcf <- function(
  vcf_dir, prefix, pattern, Genotypes = list(wt = "wildtype", mt = "mutant"), 
  min_DP = 5, min_QUAL = 5, only_mutant = FALSE) {
  
  # Find matching files
  vcf_list <- list.files(path = vcf_dir, pattern = pattern, full.names = TRUE)
  if (length(vcf_list) == 0) stop("No vcf files found.")
  
  # decide whether to load only mutant or both genotypes
  selected_genotypes <- if (only_mutant) list(mt = Genotypes[["mt"]]) else Genotypes
  geno_data <- list()
  
  # Loop through genotypes
  for (genotype in names(selected_genotypes)) {
    file_pattern <- paste0("^", prefix, "_", selected_genotypes[[genotype]])
    sample_name <- paste(prefix, selected_genotypes[[genotype]], sep = "_")
    matched_file <- vcf_list[grepl(file_pattern, basename(vcf_list))]

    if (length(matched_file) == 0) {stop(paste("File for", genotype, "not found."))}
    if (length(matched_file) > 1) {
      stop(
        paste0(
          "More than one file found for ", genotype, ".\n",
          "Please keep only ONE of the following files:\n\n",
          paste(basename(matched_file), collapse = "\n")
        )
      )
    }
    
    geno_file <- matched_file[1]
    message(paste("Reading file for", basename(geno_file), paste0("(", genotype, ").")))
    
    # Read file
    data <- data.table::fread(geno_file)
    if (is.null(data) || nrow(data) == 0) {stop(paste("File for", genotype, "is empty or invalid"))}
    message(paste(nrow(data), "rows loaded for", sample_name, paste0("(", genotype, ").")))
    
    # check and rename columns
    if (ncol(data) >= 10) {
      expected_colnames <- c("CHROM", "POS", "REF", "ALT", "QUAL", "DP", "Fref", "Rref", "Falt", "Ralt")
      data.table::setnames(data, old = colnames(data)[1:10], new = expected_colnames)
    } else {
      stop("The input file does not contain the expected number of columns (10).")
    }
    
    # Filter and Keep only SNPs (exclude indels)
    # Filter by depth and quality
    data <- data[nchar(REF) == 1 & nchar(ALT) == 1]
    data <- data[!is.na(DP) & DP > min_DP & QUAL >= min_QUAL]
    
    # Allele frequency (AF)
    data[, `:=`(Tref = Fref + Rref, Talt = Falt + Ralt)]
    data[, AF := data.table::fifelse((Tref + Talt) > 0, Talt / (Tref + Talt), NA_real_)]
    data[, nonhom := as.integer(Tref > 0 & Talt > 0)]
    
    # Sort and prefix columns
    data[, POS := as.integer(POS)]
    data.table::setorder(data, CHROM, POS)
    cols_to_prefix <- setdiff(names(data), c("CHROM", "POS"))
    data.table::setnames(data, old = cols_to_prefix, new = paste0(genotype, "_", cols_to_prefix))
    
    geno_data[[genotype]] <- data
    message(paste(nrow(data), "rows kept for", sample_name, paste0("(", genotype, ")"),"after filtering. \n"))
  }
  
  # Ensure both genotypes present if not mutant-only mode
  if (!only_mutant && (is.null(geno_data$wt) || is.null(geno_data$mt))) {
    stop("Both wild-type and mutant data are required.")
  }
  
  message(paste("Successfully created list for:", paste(names(geno_data), collapse = ", ")))
  return(geno_data)
}


#' Compute BSA metrics from bulk SNP tables
#' Computes AFD, ED, ED4, G statistics, and produces EMS + unique SNP subsets.
#' @param geno_data Output list from [read_vcf()].
#' @param prefix Prefix used for output naming.
#' @param save_results If TRUE, write tables and an RDS file.
#' @param output_dir Output directory (only used if save_results=TRUE).
#' @param only_mutant If TRUE, compute mutant-only outputs (EMS subset).
#' @return A named list of data.tables (shared/all, unique, EMS subsets).
#' @export
compute_bsa <- function(geno_data, prefix, save_results = FALSE, output_dir = "post_analysis", only_mutant = FALSE) {
  naturalsort <- function(x) {
    if (requireNamespace("stringr", quietly = TRUE)) stringr::str_sort(x, numeric = TRUE) else sort(x)
  }
  
  # function order
  sort_variants <- function(dt) {
    dt <- data.table::as.data.table(dt)
    dt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    dt[order(CHROM, POS)]
  }
  
  # ems logical helper
  is_ems <- function(ref, alt) {
    (ref == "G" & alt == "A") | (ref == "C" & alt == "T")
  }
  
  # filter and Extract EMS SNPs from mutant
  get_ems <- function(data, ref_col, alt_col) {
    ems_variants <- data[ is_ems(get(ref_col), get(alt_col)) ]
    ems_variants <- sort_variants(ems_variants)
    return(ems_variants)
  }
  
  # anti-join and identity unique SNPs
  id_unique_snps <- function(data1, data2) {
    unique_snps <- data.table::as.data.table(dplyr::anti_join(data1, data2, by = c("CHROM", "POS")))
    unique_snps <- sort_variants(unique_snps)
    return(unique_snps)
  }
  
  # Process mutant-only mode
  if (only_mutant) {
    if (is.null(geno_data$mt)) stop("Mutant data is required.")
    mt_data <- sort_variants(geno_data$mt)
    ant_mt_ems <- get_ems(mt_data, "mt_REF", "mt_ALT")
    result <- list(ant_mt = mt_data, ant_mt_ems = ant_mt_ems)
    
    if (save_results) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      saveRDS(result, file.path(output_dir, paste0(prefix, "_mutant_only_results.rds")))
      data.table::fwrite(mt_data, file.path(output_dir, paste0(prefix, "_ant_mt.csv")))
      data.table::fwrite(ant_mt_ems, file.path(output_dir, paste0(prefix, "_ant_mt_ems.csv")))
      message("Mutant-only results saved to ", output_dir)
    }
    return(result)
  } else {
    
    # === WT vs MT analysis ===
    if (is.null(geno_data$wt) || is.null(geno_data$mt)) {stop("Both wild-type and mutant datasets are required.")}
    wt_data <- data.table::as.data.table(geno_data$wt)
    mt_data <- data.table::as.data.table(geno_data$mt)
    
    # Merge wild-type and mutant by CHROM and POS
    message("Merging wild-type and mutant datasets...")
    wt_mt <- merge(wt_data, mt_data, by = c("CHROM", "POS"))
    wt_mt <- sort_variants(wt_mt)
    message("Number of shared SNPs (wt_mt): ", nrow(wt_mt))
    
    # Statistical metrics : Allele frequency diff, ED, G-statistics, p-values
    wt_mt[, AFD := abs(wt_AF - mt_AF)]
    
    # Calculate Euclidean Distance (ED) and its fourth power (ED4)
    wt_mt[, ED := sqrt(2) * AFD]
    wt_mt[, ED4 := ED^4]
    
    # Calculate G-statistics for each SNP
    compute_G <- function(wt_Tref, wt_Talt, mt_Tref, mt_Talt) {
      # Compute allele counts
      mt_ref <- mt_Tref # mt_Fref + mt_Rref   (reference allele in mutant bulk)
      mt_alt <- mt_Talt # mt_Falt + mt_Ralt   (alternate allele in mutant bulk)
      wt_ref <- wt_Tref # wt_Fref + wt_Rref   (reference allele in wildtype bulk)
      wt_alt <- wt_Talt # wt_Falt + wt_Ralt   (alternate allele in wildtype bulk)
      
      total <- mt_ref + mt_alt + wt_ref + wt_alt
      if (total == 0) return(NA_real_)
      
      row_ref <- mt_ref + wt_ref
      row_alt <- mt_alt + wt_alt
      col_mt <- mt_ref + mt_alt
      col_wt <- wt_ref + wt_alt
      
      # Expected counts under independence
      e_mt_ref <- col_mt * row_ref / total
      e_mt_alt <- col_mt * row_alt / total
      e_wt_ref <- col_wt * row_ref / total
      e_wt_alt <- col_wt * row_alt / total
      
      obs <- c(mt_ref, mt_alt, wt_ref, wt_alt)
      exp <- c(e_mt_ref, e_mt_alt, e_wt_ref, e_wt_alt)
      
      # Avoid division by zero or log(0)
      obs[obs == 0] <- 1e-10
      exp[exp == 0] <- 1e-10
      
      return(2 * sum(obs * log(obs / exp), na.rm = TRUE))
    }
    
    wt_mt[, G := mapply(compute_G, wt_Tref, wt_Talt, mt_Tref, mt_Talt)]
    
    
    # Unique SNPs and EMS filtering
    ant_wt <- id_unique_snps(wt_data, mt_data)
    ant_mt <- id_unique_snps(mt_data, wt_data)
    ant_wt_ems <- get_ems(ant_wt, "wt_REF", "wt_ALT")
    ant_mt_ems <- get_ems(ant_mt, "mt_REF", "mt_ALT")
    wt_mt_ems <- wt_mt[ is_ems(wt_REF, wt_ALT) & is_ems(mt_REF, mt_ALT) ]
    wt_mt_ems <- sort_variants(wt_mt_ems)
    
    result <- list(
      wt_mt_all = wt_mt, ant_mt_all = ant_mt, ant_wt_all = ant_wt, 
      wt_mt_ems = wt_mt_ems, ant_mt_ems = ant_mt_ems, ant_wt_ems = ant_wt_ems
      )
    
    # === Save output ===
    if (save_results) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      saveRDS(result, file.path(output_dir, paste0(prefix, "_results.rds")))
      for (name in names(result)) {
        data.table::fwrite(result[[name]], file.path(output_dir, paste0(prefix, "_", name, ".csv")))
      }
      message("All results to ", output_dir)
    }
    return(result)
  }
}