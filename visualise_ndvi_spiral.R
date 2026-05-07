# ==============================================================================
# NDVI Seasonality Spiral — Wild Ennerdale (demo site)
# Manual polar→Cartesian with arc interpolation for smooth curves
# ==============================================================================

library(terra)
library(tidyverse)
library(cowplot)

# ══════════════════════════════════════════════════════════
#  FUNCTIONS
# ══════════════════════════════════════════════════════════

load_ndvi_long <- function(input_file) {
  r <- rast(input_file)
  NAflag(r) <- -9999
  dates <- as.Date(gsub("^NDVI_", "", names(r)), format = "%Y%m%d")
  vals <- values(r)
  colnames(vals) <- as.character(dates)
  valid_px <- rowSums(!is.na(vals)) > 0
  vals <- vals[valid_px, , drop = FALSE]
  df <- as.data.frame(vals)
  df$pixel_id <- seq_len(nrow(df))
  df %>%
    pivot_longer(cols = -pixel_id, names_to = "date_str", values_to = "ndvi") %>%
    mutate(
      date = as.Date(date_str),
      year = lubridate::year(date),
      doy  = lubridate::yday(date)
    ) %>%
    filter(!is.na(ndvi), ndvi >= 0, ndvi <= 1) %>%
    select(pixel_id, date, year, doy, ndvi) %>%
    arrange(pixel_id, date)
}

polar_to_cartesian <- function(df, jitter_days = 0) {
  # Convert doy + ndvi to Cartesian (x, y) for manual polar plot
  # Jan 1 at bottom (6 o'clock), clockwise
  # jitter_days: uniform noise ± this many days, to spread shared composite dates
  df %>%
    mutate(
      doy_j = doy + runif(n(), -jitter_days, jitter_days),
      angle = -pi/2 - 2 * pi * doy_j / 365.25,
      cx = ndvi * cos(angle),
      cy = ndvi * sin(angle)
    )
}

interpolate_arcs <- function(df, n_interp = 5, max_gap_days = 50) {
  # For each pixel, interpolate between consecutive observations
  # so that the path follows the circle instead of cutting chords.
  # Each pair of consecutive obs gets n_interp intermediate points.
  # Interpolation is done in polar (doy, ndvi) space, then converted to Cartesian.
  #
  # For Dec→Jan transitions (where doy decreases), the interpolation
  # wraps through 365 naturally.
  #

  # Gaps longer than max_gap_days insert an NA break instead of interpolating,
  # preventing spurious arcs through low-NDVI territory during winter data gaps.

  df %>%
    group_by(pixel_id) %>%
    group_modify(function(grp, key) {
      n <- nrow(grp)
      if (n < 2) return(tibble())

      out <- vector("list", n - 1)
      for (i in seq_len(n - 1)) {
        gap_days <- as.numeric(grp$date[i + 1] - grp$date[i])

        # If the gap is too large, insert a break instead of interpolating
        if (gap_days > max_gap_days) {
          out[[i]] <- tibble(
            doy  = c(grp$doy[i], NA),
            ndvi = c(grp$ndvi[i], NA),
            year = c(grp$year[i], NA)
          )
          next
        }

        doy1  <- grp$doy[i]
        doy2  <- grp$doy[i + 1]
        ndvi1 <- grp$ndvi[i]
        ndvi2 <- grp$ndvi[i + 1]
        yr1   <- grp$year[i]
        yr2   <- grp$year[i + 1]

        # Handle wrap: if doy decreases, we're crossing Dec→Jan
        if (doy2 < doy1) {
          doy2_adj <- doy2 + 365.25
        } else {
          doy2_adj <- doy2
        }

        t_seq <- seq(0, 1, length.out = n_interp + 2)
        doy_interp  <- doy1 + t_seq * (doy2_adj - doy1)
        ndvi_interp <- ndvi1 + t_seq * (ndvi2 - ndvi1)
        year_interp <- yr1 + t_seq * (yr2 - yr1)

        # Wrap doy back into [0, 365.25]
        doy_interp <- doy_interp %% 365.25

        out[[i]] <- tibble(
          doy  = doy_interp,
          ndvi = ndvi_interp,
          year = year_interp
        )
      }
      bind_rows(out)
    }) %>%
    ungroup() %>%
    mutate(
      angle = -pi/2 - 2 * pi * doy / 365.25,
      cx = ndvi * cos(angle),
      cy = ndvi * sin(angle)
    )
}

make_grid <- function() {
  ring_vals <- c(0.2, 0.4, 0.6, 0.8)
  theta_seq <- seq(0, 2 * pi, length.out = 200)
  rings <- map_dfr(ring_vals, function(rv) {
    tibble(x = rv * cos(theta_seq), y = rv * sin(theta_seq), r = rv)
  })

  month_doys <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
  month_labs <- c("Jan","Feb","Mar","Apr","May","Jun",
                  "Jul","Aug","Sep","Oct","Nov","Dec")
  spokes <- tibble(doy = month_doys, label = month_labs) %>%
    mutate(
      angle = -pi/2 - 2 * pi * doy / 365.25,
      x0 = 0, y0 = 0,
      x1 = 1.02 * cos(angle), y1 = 1.02 * sin(angle),
      xl = 1.12 * cos(angle), yl = 1.12 * sin(angle)
    )

  ring_labels <- tibble(r = ring_vals, x = 0.02, y = ring_vals,
                         label = as.character(ring_vals))

  list(rings = rings, spokes = spokes, ring_labels = ring_labels)
}

build_spiral_plot <- function(path_df, grid, yr_range,
                               alpha = 0.012, linewidth = 0.2,
                               point_size = 0.15, use_points = FALSE,
                               title_suffix = "") {
  p <- ggplot() +
    geom_path(data = grid$rings, aes(x = x, y = y, group = r),
              colour = "black", linewidth = 0.25, alpha = 0.3) +
    geom_segment(data = grid$spokes,
                 aes(x = x0, y = y0, xend = x1, yend = y1),
                 colour = "black", linewidth = 0.25, alpha = 0.3)

  if (use_points) {
    p <- p + geom_point(data = path_df,
                        aes(x = cx, y = cy, colour = year),
                        alpha = alpha, size = point_size, shape = 16)
  } else {
    p <- p + geom_path(data = path_df,
                       aes(x = cx, y = cy, colour = year, group = pixel_id),
                       alpha = alpha, linewidth = linewidth)
  }

  p <- p +
    geom_text(data = grid$spokes,
              aes(x = xl, y = yl, label = label),
              size = 2.5, colour = "black", family = "Helvetica",
              fontface = "bold") +
    geom_text(data = grid$ring_labels,
              aes(x = x, y = y, label = label),
              size = 1.8, colour = "black", family = "Helvetica",
              hjust = 1, vjust = -0.3) +
    scale_colour_viridis_c(
      name = NULL, option = "inferno", begin = 0.15, end = 0.9,
      breaks = pretty(yr_range[1]:yr_range[2], n = 5),
      guide = guide_colourbar(barwidth = unit(0.3, "cm"),
                               barheight = unit(4, "cm"), ticks = FALSE)
    ) +
    coord_fixed(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
    labs(
      title = title_suffix,
      subtitle = NULL
    ) +
    theme_void(base_size = 9, base_family = "Helvetica") +
    theme(
      plot.title      = element_text(size = 10, face = "bold", hjust = 0.5,
                                      margin = margin(b = 2)),
      plot.subtitle   = element_text(size = 7.5, hjust = 0.5, colour = "grey40",
                                      margin = margin(b = 8)),
      legend.position.inside = c(0.90, 0.15),
      legend.text     = element_text(size = 6),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin     = margin(10, 10, 10, 10)
    )
  p
}

# ══════════════════════════════════════════════════════════
#  ANNUAL METRIC FUNCTIONS
# ══════════════════════════════════════════════════════════

load_annual_metrics <- function(mk_file, metric, start_year = 2003,
                                 end_year = 2022) {
  # Extract a single annual metric (INDVI, minNDVI, maxNDVI) from the stacked
  # MK/Annual GeoTIFF and return a long-format tibble.
  r <- rast(mk_file)
  NAflag(r) <- -9999
  years <- start_year:end_year
  band_names <- paste0(metric, "_", years)
  r_sub <- r[[band_names]]
  vals <- values(r_sub)
  colnames(vals) <- as.character(years)
  valid_px <- rowSums(!is.na(vals)) > 0
  vals <- vals[valid_px, , drop = FALSE]
  df <- as.data.frame(vals)
  df$pixel_id <- seq_len(nrow(df))
  df %>%
    pivot_longer(-pixel_id, names_to = "year", values_to = "value") %>%
    mutate(year = as.integer(year), metric = metric) %>%
    filter(!is.na(value))
}

build_trend_panel <- function(df, metric_label, y_label,
                               yr_range, colour_hex = "grey30") {
  # PNAS-style trend panel: per-pixel semi-transparent lines + loess mean
  ggplot(df, aes(x = year, y = value)) +
    geom_line(aes(group = pixel_id), alpha = 0.03, linewidth = 0.3,
              colour = colour_hex) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                colour = "black", fill = "grey70",
                linewidth = 0.6, alpha = 0.3) +
    scale_x_continuous(breaks = seq(yr_range[1], yr_range[2], by = 4),
                       limits = yr_range, expand = c(0.01, 0)) +
    labs(x = NULL, y = y_label, title = metric_label) +
    theme_classic(base_size = 8, base_family = "Helvetica") +
    theme(
      plot.title      = element_text(size = 8, face = "bold", hjust = 0,
                                      margin = margin(b = 2)),
      axis.title.y    = element_text(size = 7),
      axis.text       = element_text(size = 6, colour = "black"),
      axis.line       = element_line(linewidth = 0.3),
      axis.ticks      = element_line(linewidth = 0.3),
      plot.margin     = margin(4, 8, 4, 4)
    )
}

# ══════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════

generate_spiral_composite <- function(ts_file, mk_file, output_dir,
                                      site_name = "Wild Ennerdale") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ── Load time series data ────────────────────────────────
df_long  <- load_ndvi_long(ts_file)
yr_range <- range(df_long$year)
grid     <- make_grid()

# ── Panel A: Spiral (all pixels, jittered points) ────────
all_cart <- polar_to_cartesian(df_long, jitter_days = 8)
p_spiral <- build_spiral_plot(all_cart, grid, yr_range,
                               alpha = 0.08, use_points = TRUE,
                               point_size = 0.4)

# ── Load annual metrics ──────────────────────────────────
df_indvi  <- load_annual_metrics(mk_file, "INDVI")
df_min    <- load_annual_metrics(mk_file, "minNDVI")
df_max    <- load_annual_metrics(mk_file, "maxNDVI")

# ── Panels B–D: Annual trend plots ───────────────────────
  p_indvi <- build_trend_panel(df_indvi, "B",
                                "iNDVI", yr_range, "#2a6e3f")
  p_min   <- build_trend_panel(df_min, "C",
                                "min NDVI", yr_range, "#8c510a")
  p_max   <- build_trend_panel(df_max, "D",
                                "max NDVI", yr_range, "#01665e")

# ── Composite figure ─────────────────────────────────────
# Left: spiral (full height)  |  Right: 3 stacked trend panels
right_col <- plot_grid(p_indvi, p_min, p_max,
                        ncol = 1, align = "v", axis = "lr")

p_spiral_labelled <- p_spiral +
  labs(title = "A") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0,
                                    margin = margin(b = 2)))

composite <- plot_grid(p_spiral_labelled, right_col,
                        ncol = 2, rel_widths = c(1, 0.85))

out_composite <- file.path(output_dir, paste0("NDVI_spiral_composite_", gsub(" ", "_", site_name), ".png"))
ggsave(out_composite, composite,
       width = 10, height = 5.5, dpi = 300, bg = "white")
message(paste("\u2705 Composite figure saved to", out_composite))

# ── Also save standalone spiral ──────────────────────────
out_standalone <- file.path(output_dir, paste0("NDVI_spiral_", gsub(" ", "_", site_name), ".png"))
ggsave(out_standalone, p_spiral, width = 6, height = 6, dpi = 300, bg = "white")
message(paste("\u2705 Standalone spiral saved to", out_standalone))

# ── Single pixel diagnostic ─────────────────────────────
best_px <- df_long %>% count(pixel_id) %>% slice_max(n, n = 1) %>%
  pull(pixel_id) %>% `[`(1)
  single_path <- df_long %>% filter(pixel_id == best_px) %>%
    interpolate_arcs(n_interp = 8, max_gap_days = 50)
  p_single <- build_spiral_plot(single_path, grid, yr_range,
                                 alpha = 0.8, linewidth = 0.5,
                                 title_suffix = " (single pixel)")
  out_single <- file.path(output_dir, paste0("NDVI_spiral_", gsub(" ", "_", site_name), "_single_pixel.png"))
  ggsave(out_single, p_single, width = 6, height = 6, dpi = 300, bg = "white")
  message(paste("\u2705 Single-pixel spiral saved to", out_single))
}

# ══════════════════════════════════════════════════════════
#  EXECUTION
# ══════════════════════════════════════════════════════════

if (sys.nframe() == 0) {
  generate_spiral_composite(
    ts_file = "Outputs/GEE_MK_250m/NDVI_TimeSeries_250m_Wild_Ennerdale.tif",
    mk_file = "Outputs/GEE_MK_250m/MK_Stacked_250m_Wild_Ennerdale.tif",
    output_dir = "Outputs/plots",
    site_name = "Wild Ennerdale"
  )
}

