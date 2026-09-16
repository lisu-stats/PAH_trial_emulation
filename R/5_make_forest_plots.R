########################################################################
## Forest Plots for MSM Results (Main, Treat×Year, GRIPHON eligibility)
########################################################################

## ===== 0) Libraries ==========================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

## ===== 1) Output Settings ==============================================
IMAGE_FORMAT <- "tiff"  # Options: "tiff", "jpeg", "jpg", "png"
IMAGE_DPI    <- 300     # High resolution for Word documents and publication

## ===== 2) Plot helpers =======================================================

# 2A. Forest plot for risk (never-initiate) by year
plot_forest_risk <- function(df,
                             x_year = "horizon",
                             y_risk = "risk_never",
                             lo = "risk_never_lo",
                             hi = "risk_never_hi",
                             title = NULL,
                             subtitle = NULL,
                             xlab = "Risk (never-initiate)",
                             limits = NULL,
                             accuracy = 1) {
  
  df <- as.data.frame(df)
  req <- c(x_year, y_risk, lo, hi)
  stopifnot(all(req %in% names(df)))
  
  year_levels <- rev(unique(df[[x_year]]))
  
  df %>%
    mutate(.Year = factor(.data[[x_year]], levels = year_levels)) %>%
    ggplot(aes(x = .data[[y_risk]], y = .Year)) +
    geom_errorbar(
      aes(xmin = .data[[lo]], xmax = .data[[hi]]),
      orientation = "y",
      width = 0.2,
      linewidth = 0.7
    ) +
    geom_point(size = 2.7) +
    scale_x_continuous(
      labels = percent_format(accuracy = accuracy), 
      limits = limits,
      expand = expansion(mult = c(0.05, 0.15))
    ) +
    labs(x = xlab, y = "Year", title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
    )
}

# 2B. Forest plot for risk differences (initiate − never) in Percentages
plot_forest_rd <- function(df,
                           x_year = "horizon",
                           y_rd = "RD",
                           lo = "RD_lo",
                           hi = "RD_hi",
                           title = NULL,
                           subtitle = NULL,
                           xlab = "Risk difference (initiate − never)",
                           accuracy = 1) {
  
  df <- as.data.frame(df)
  req <- c(x_year, y_rd, lo, hi)
  stopifnot(all(req %in% names(df)))
  
  year_levels <- rev(unique(df[[x_year]]))
  
  df %>%
    mutate(.Year = factor(.data[[x_year]], levels = year_levels)) %>%
    ggplot(aes(x = .data[[y_rd]], y = .Year)) +
    geom_vline(xintercept = 0, linetype = 3, color = "grey40") +
    geom_errorbar(
      aes(xmin = .data[[lo]], xmax = .data[[hi]]),
      orientation = "y",
      width = 0.2,
      linewidth = 0.7
    ) +
    geom_point(size = 2.7) +
    scale_x_continuous(
      labels = percent_format(accuracy = accuracy),
      expand = expansion(mult = c(0.08, 0.15))
    ) +
    labs(x = xlab, y = "Year", title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
    )
}

# 2C. Forest plot for MSM coefficients (logit scale)
plot_forest_coef <- function(df,
                             coef_col = "coefficient",
                             est = "estimate",
                             lo = "CI_lo",
                             hi = "CI_hi",
                             title = NULL,
                             subtitle = NULL,
                             xlab = "MSM coefficient (logit scale)") {
  
  df <- as.data.frame(df)
  req <- c(coef_col, est, lo, hi)
  stopifnot(all(req %in% names(df)))
  
  level_order <- rev(df[[coef_col]])
  
  df %>%
    mutate(.Param = factor(.data[[coef_col]], levels = level_order)) %>%
    ggplot(aes(x = .data[[est]], y = .Param)) +
    geom_vline(xintercept = 0, linetype = 3, color = "grey40") +
    geom_errorbar(
      aes(xmin = .data[[lo]], xmax = .data[[hi]]),
      orientation = "y",
      width = 0.2,
      linewidth = 0.7
    ) +
    geom_point(size = 2.7) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.15))) +
    labs(x = xlab, y = NULL, title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
    )
}

# 2D. Faceted RD plot for strata in Percentages
plot_forest_rd_facet <- function(df,
                                 x_year = "horizon",
                                 y_rd = "RD",
                                 lo = "RD_lo",
                                 hi = "RD_hi",
                                 title = NULL,
                                 subtitle = NULL,
                                 accuracy = 1) {
  
  df <- as.data.frame(df)
  req <- c(x_year, y_rd, lo, hi, "stratum")
  stopifnot(all(req %in% names(df)))
  
  year_levels <- rev(unique(df[[x_year]]))
  
  df %>%
    mutate(.Year = factor(.data[[x_year]], levels = year_levels)) %>%
    ggplot(aes(x = .data[[y_rd]], y = .Year)) +
    geom_vline(xintercept = 0, linetype = 3, color = "grey40") +
    geom_errorbar(
      aes(xmin = .data[[lo]], xmax = .data[[hi]]),
      orientation = "y",
      width = 0.2,
      linewidth = 0.7
    ) +
    geom_point(size = 2.7) +
    facet_wrap(~ stratum, ncol = 1, scales = "fixed") +
    scale_x_continuous(
      labels = percent_format(accuracy = accuracy),
      expand = expansion(mult = c(0.08, 0.15))
    ) +
    labs(x = "Risk difference (initiate − never)", y = "Year",
         title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
    )
}

# 2E. Flexible image saver (TIFF / JPEG / PNG)
save_plot <- function(plot, base_filename, width = 6, height = 3.6, 
                      format = IMAGE_FORMAT, dpi = IMAGE_DPI) {
  
  ext <- paste0(".", gsub("^\\.", "", format))
  filename <- paste0(base_filename, ext)
  
  message("Saving: ", filename, " (DPI = ", dpi, ")")
  
  if (tolower(format) %in% c("tif", "tiff")) {
    ggsave(filename = filename, plot = plot, width = width, height = height, 
           dpi = dpi, compression = "lzw")
  } else {
    ggsave(filename = filename, plot = plot, width = width, height = height, 
           dpi = dpi)
  }
}

## ===== 3) Small utilities ====================================================

coerce_cols_rd <- function(df) {
  df <- as.data.frame(df)
  if (ncol(df) >= 9) {
    df <- df[, 1:9]
    colnames(df) <- c("horizon", "risk_never", "risk_never_SE", 
                      "risk_never_lo", "risk_never_hi", "RD", 
                      "RD_SE", "RD_lo", "RD_hi")
  }
  needed <- c("horizon", "risk_never", "risk_never_lo", "risk_never_hi", "RD", "RD_lo", "RD_hi")
  miss <- setdiff(needed, names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  df[, !is.na(names(df)) & names(df) != "", drop = FALSE]
}

coerce_cols_coef <- function(df) {
  df <- as.data.frame(df)
  if (ncol(df) >= 4) {
    df <- df[, 1:4]
    colnames(df) <- c("coefficient", "estimate", "CI_lo", "CI_hi")
  }
  needed <- c("coefficient", "estimate", "CI_lo", "CI_hi")
  miss <- setdiff(needed, names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  df[, !is.na(names(df)) & names(df) != "", drop = FALSE]
}

## ===== 4) Theme Settings ====================================================
theme_set(theme_minimal(base_size = 12))
update_geom_defaults("errorbar", list(linewidth = 0.7))
update_geom_defaults("point", list(size = 2.7))

## ===== 5) Map Object Names ==================================================
if (exists("RD_summary_main_g20")) {
  MAIN_RD_G20 <- coerce_cols_rd(RD_summary_main_g20)
}

if (exists("RD_summary_trtyear_g20")) {
  TBY_RD_G20 <- coerce_cols_rd(RD_summary_trtyear_g20)
}

if (exists("RD_summary_GRIPHON_g20")) {
  GRIPHON_RD_G0_G20 <- coerce_cols_rd(RD_summary_GRIPHON_g20[, 1:9])
  GRIPHON_RD_G1_G20 <- coerce_cols_rd(RD_summary_GRIPHON_g20[, c(1, 10:17)])
}

if (exists("RD_summary_prevalent_g20")) {
  prevalent_RD_G0_G20 <- coerce_cols_rd(RD_summary_prevalent_g20[, 1:9])
  prevalent_RD_G1_G20 <- coerce_cols_rd(RD_summary_prevalent_g20[, c(1, 10:17)])
}

## ===== 6) Plotting Pipeline =================================================

if (exists("MAIN_RD_G20")) {
  p <- plot_forest_risk(MAIN_RD_G20)
  save_plot(p, "fig_main_risk_g20")
  p <- plot_forest_rd(MAIN_RD_G20)
  save_plot(p, "fig_main_rd_g20")
}

if (exists("TBY_RD_G20")) {
  p <- plot_forest_risk(TBY_RD_G20)
  save_plot(p, "fig_tby_risk_g20")
  p <- plot_forest_rd(TBY_RD_G20)
  save_plot(p, "fig_tby_rd_g20")
}

if (exists("GRIPHON_RD_G0_G20")) {
  p <- plot_forest_risk(GRIPHON_RD_G0_G20)
  save_plot(p, "fig_griphon_g0_risk_g20")
  p <- plot_forest_rd(GRIPHON_RD_G0_G20)
  save_plot(p, "fig_griphon_g0_rd_g20")
}

if (exists("GRIPHON_RD_G1_G20")) {
  p <- plot_forest_risk(GRIPHON_RD_G1_G20)
  save_plot(p, "fig_griphon_g1_risk_g20")
  p <- plot_forest_rd(GRIPHON_RD_G1_G20)
  save_plot(p, "fig_griphon_g1_rd_g20")
}

if (exists("GRIPHON_RD_G0_G20") && exists("GRIPHON_RD_G1_G20")) {
  df_fac_rd_g20 <- dplyr::bind_rows(
    GRIPHON_RD_G0_G20 %>% dplyr::mutate(stratum = "GRIPHON = No"),
    GRIPHON_RD_G1_G20 %>% dplyr::mutate(stratum = "GRIPHON = Yes")
  )
  p <- plot_forest_rd_facet(df_fac_rd_g20)
  save_plot(p, "fig_griphon_rd_faceted_g20", height = 5.0)
  
  df_fac_risk_g20 <- dplyr::bind_rows(
    GRIPHON_RD_G0_G20 %>% dplyr::mutate(stratum = "GRIPHON = No"),
    GRIPHON_RD_G1_G20 %>% dplyr::mutate(stratum = "GRIPHON = Yes")
  )
  p <- plot_forest_risk(df_fac_risk_g20) + facet_wrap(~ stratum, ncol = 1)
  save_plot(p, "fig_griphon_risk_faceted_g20", height = 5.0)
}

if (exists("prevalent_RD_G0_G20") && exists("prevalent_RD_G1_G20")) {
  df_fac_rd_g20 <- dplyr::bind_rows(
    prevalent_RD_G0_G20 %>% dplyr::mutate(stratum = "prevalent = No"),
    prevalent_RD_G1_G20 %>% dplyr::mutate(stratum = "prevalent = Yes")
  )
  p <- plot_forest_rd_facet(df_fac_rd_g20)
  save_plot(p, "fig_prevalent_rd_faceted_g20", height = 5.0)
  
  df_fac_risk_g20 <- dplyr::bind_rows(
    prevalent_RD_G0_G20 %>% dplyr::mutate(stratum = "prevalent = No"),
    prevalent_RD_G1_G20 %>% dplyr::mutate(stratum = "prevalent = Yes")
  )
  p <- plot_forest_risk(df_fac_risk_g20) + facet_wrap(~ stratum, ncol = 1)
  save_plot(p, "fig_prevalent_risk_faceted_g20", height = 5.0)
}


if (exists("RD_summary_RPH_g20")) {
  RPH_RD_G0_G20 <- coerce_cols_rd(RD_summary_RPH_g20[, 1:9])
  RPH_RD_G1_G20 <- coerce_cols_rd(RD_summary_RPH_g20[, c(1, 10:17)])
}


if (exists("RPH_RD_G0_G20") && exists("RPH_RD_G1_G20")) {
  df_fac_rd_g20 <- dplyr::bind_rows(
    RPH_RD_G0_G20 %>% dplyr::mutate(stratum = "RPH = No"),
    RPH_RD_G1_G20 %>% dplyr::mutate(stratum = "RPH = Yes")
  )
  p <- plot_forest_rd_facet(df_fac_rd_g20)
  save_plot(p, "fig_RPH_rd_faceted_g20", height = 5.0)
  
  df_fac_risk_g20 <- dplyr::bind_rows(
    RPH_RD_G0_G20 %>% dplyr::mutate(stratum = "RPH = No"),
    RPH_RD_G1_G20 %>% dplyr::mutate(stratum = "RPH = Yes")
  )
  p <- plot_forest_risk(df_fac_risk_g20) + facet_wrap(~ stratum, ncol = 1)
  save_plot(p, "fig_RPH_risk_faceted_g20", height = 5.0)
}
