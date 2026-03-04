############################################################
# R/03_plot.R
# Plotting functions:
#   - bsa_theme()
#   - plot_bsa()
############################################################

bsa_theme <- function(base_size = 60, legendit = "none") {
  require(grid)
  theme_linedraw(base_size = base_size) + 
    theme(
      plot.title = element_text(size = base_size, face = "bold", hjust = 0.5, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size * 0.7, face = "bold", color = "black"),
      axis.text.y = element_text(size = base_size, face = "bold", colour = "black", hjust = 1),
      axis.title.x = element_text(size = base_size, face = "bold", color = "black"),
      axis.title.y = element_text(size = base_size, face = "bold", angle = 90, color = "black"),
      axis.line = element_line(colour = "black", linewidth = 3),
      axis.ticks = element_line(colour = "black", linewidth = 4),
      axis.ticks.length = unit(0.5, "cm"),
      panel.border = element_rect(linewidth = 7, fill = NA),
      panel.spacing.x = unit(2, "lines"),
      panel.spacing.y = unit(2, "lines"),
      panel.grid.minor = element_line(colour = "grey90", linewidth = 2),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.0),
      strip.background = element_rect(fill = "grey80", colour = "black", linewidth = 7),
      strip.text = element_text(colour = "black", face = "bold", size = base_size),
      strip.text.x = element_text(colour = "black", face = "bold", size = base_size, margin = margin(0.6, 0.0, 0.6, 0.0, "cm")),
      strip.text.y = element_text(colour = "black", face = "bold", size = base_size, margin = margin(0.0, 0.6, 0.0, 0.6, "cm")), 
      legend.position = legendit
    )
}


# Plot BSA-Seq Data by Chromosome
plot_bsa <- function(
    data, prefix, column, y_title, plot_title, file_suffix, 
    ylim = NULL, is_locfit = FALSE, is_rollmedian = TRUE, bwidth = 1000000,
    is_histogram = FALSE, threshold = NULL, rollmedian = 501,
    hwidth, hheight, width, height, dpi, device, plots_dir, 
    nn_prop, point_size = 4, line_size = 4, alpha_size = 1, facet_column = 5, 
    plot_style = c("wrap", "grid"), remove_x_text = FALSE, color_panel = c("blue", "red")) {
  
  # Ensure plots directory exists
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # Helper mini-function : Capitalize the first letter 
  capitalize_first <- function(text) {
    substr(text, 1, 1) <- toupper(substr(text, 1, 1))
    text
  }
  
  # Capitalize inbred name from prefix
  inbred <- capitalize_first(prefix)
  
  # Add PositionMb column for better x-axis readability
  data <- data %>% dplyr::mutate(PositionMb = POS / 1e6)
  bwidth <- bwidth / 1e6
  
  # Define color panel for chromosomes
  color_panel <- rep(color_panel, length.out = length(unique(data$CHROM)))
  
  # If smooth is requested then rollmedian is irrelevant
  if (is_locfit) is_rollmedian <- FALSE
  
  # plot_style
  plot_style <- match.arg(plot_style)
  is_wrap <- plot_style == "wrap"
  
  # Select faceting style
  facet_layer <- if (is_wrap) {
    facet_wrap(~ CHROM, ncol = facet_column, scales = "free_x")
  } else {
    facet_grid(. ~ CHROM, scales = "free_x", space = "free_x")
  }
  
  # Shared labels + theme bits
  title_txt <- paste0("Aligned to ", inbred, " :  ", plot_title)
  xlab_txt  <- "Chromosome Position (Mb) \n"
  ylab_txt  <- paste0("\n", y_title, "\n")
  
  axis_x_theme <- theme(
    axis.text.x  = if (remove_x_text) element_blank() else element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = -1)),
    axis.ticks.x = if (remove_x_text) element_blank() else element_line(), axis.title.x = element_text(color = "black")
  )
  
  # Base plot (shared scaffolding)
  plot <- ggplot(data) + facet_layer + labs(title = title_txt, x = xlab_txt, y = ylab_txt) + bsa_theme() + axis_x_theme
  
  # Build plot
  if (!is_histogram) {
    plot <- plot +
      aes(x = PositionMb, y = !!rlang::sym(column), color = CHROM) +  scale_color_manual(values = color_panel) + guides(color = "none")
    
    # Layers
    if (is_locfit) {
      plot <- plot + stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop), linewidth = line_size)
    } else if (is_rollmedian) {
      plot <- plot + geom_line(aes(y = zoo::rollmedian(!!rlang::sym(column), k = rollmedian, na.pad = TRUE)), linewidth = line_size)
    } else {
      plot <- plot + geom_point(size = point_size)
    }
    
    # Threshold (optional)
    if (!is.null(threshold)) {
      plot <- plot + geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = line_size)
    }
    
  } else {
    plot <- plot +
      aes(x = PositionMb, fill = CHROM) + geom_histogram(binwidth = bwidth, alpha = alpha_size, linewidth = 6) +
      scale_fill_manual(values = color_panel) + scale_y_continuous(labels = scales::label_number(scale = 1e-3, accuracy = 0.01, trim = TRUE))
  }
  
  # Apply fixed y-axis limits only if ylim is not NULL
  if (!is.null(ylim)) plot <- plot + coord_cartesian(ylim = ylim)
  
  # date tag for filenames (YYYYMMDD)
  date_tag <- format(Sys.Date(), "%Y%m%d")
  
  # Construct file path and save plot
  file_path   <- file.path(plots_dir, paste0(date_tag, "_", inbred, "_", file_suffix, "_", plot_style, ".", device))
  plot_width  <- if (is_wrap) hwidth else width
  plot_height <- if (is_wrap) hheight else height
  
  ggsave(filename = file_path, plot = plot, device = device, width = plot_width, height = plot_height, dpi = dpi)
  
  message(paste0("Plot saved to: ", file_path))
  invisible(plot)
}
