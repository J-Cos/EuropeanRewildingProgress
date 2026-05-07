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

interpolate_arcs <- function(df, n_interp = 5) {
  # For each pixel, interpolate between consecutive observations
  # so that the path follows the circle instead of cutting chords.
  # Each pair of consecutive obs gets n_interp intermediate points.
  # Interpolation is done in polar (doy, ndvi) space, then converted to Cartesian.
  #
  # For Dec→Jan transitions (where doy decreases), the interpolation
  # wraps through 365 naturally.

  df %>%
    group_by(pixel_id) %>%
    group_modify(function(grp, key) {
      n <- nrow(grp)
      if (n < 2) return(tibble())

      # Build interpolated points between each consecutive pair
      out <- vector("list", n - 1)
      for (i in seq_len(n - 1)) {
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
      angle = -pi/2 + 2 * pi * doy / 365.25,
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
      angle = -pi/2 + 2 * pi * doy / 365.25,
      x0 = 0, y0 = 0,
      x1 = 1.02 * cos(angle), y1 = 1.02 * sin(angle),
      xl = 1.12 * cos(angle), yl = 1.12 * sin(angle)
    )

  ring_labels <- tibble(r = ring_vals, x = -0.02, y = ring_vals,
                         label = as.character(ring_vals))

  list(rings = rings, spokes = spokes, ring_labels = ring_labels)
}

build_spiral_plot <- function(path_df, grid, yr_range,
                               alpha = 0.012, linewidth = 0.2,
                               title_suffix = "") {
  ggplot() +
    geom_path(data = grid$rings, aes(x = x, y = y, group = r),
              colour = "grey88", linewidth = 0.3) +
    geom_segment(data = grid$spokes,
                 aes(x = x0, y = y0, xend = x1, yend = y1),
                 colour = "grey88", linewidth = 0.3) +
    geom_path(data = path_df,
              aes(x = cx, y = cy, colour = year, group = pixel_id),
              alpha = alpha, linewidth = linewidth) +
    geom_text(data = grid$spokes,
              aes(x = xl, y = yl, label = label),
              size = 2.5, colour = "grey30", family = "Helvetica") +
    geom_text(data = grid$ring_labels,
              aes(x = x, y = y, label = label),
              size = 1.8, colour = "grey50", family = "Helvetica",
              hjust = 1, vjust = -0.3) +
    scale_colour_viridis_c(
      name = NULL, option = "inferno", begin = 0.15, end = 0.9,
      breaks = pretty(yr_range[1]:yr_range[2], n = 5),
      guide = guide_colourbar(barwidth = unit(0.3, "cm"),
                               barheight = unit(4, "cm"), ticks = FALSE)
    ) +
    coord_fixed(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
    labs(
      title = paste0("Wild Ennerdale \u2014 NDVI Seasonality", title_suffix),
      subtitle = paste0("(", yr_range[1], "\u2013", yr_range[2], ")")
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
}

# ══════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════

input_file <- "Outputs/GEE_MK_250m/NDVI_TimeSeries_250m_Wild_Ennerdale.tif"
output_dir <- "Outputs/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df_long  <- load_ndvi_long(input_file)
yr_range <- range(df_long$year)
grid     <- make_grid()

# ── Single pixel diagnostic ─────────────────────────────
best_px <- df_long %>% count(pixel_id) %>% slice_max(n, n = 1) %>%
  pull(pixel_id) %>% `[`(1)
cat("Diagnostic pixel:", best_px, "\n")

single_path <- df_long %>% filter(pixel_id == best_px) %>%
  interpolate_arcs(n_interp = 8)
p_single <- build_spiral_plot(single_path, grid, yr_range,
                               alpha = 0.8, linewidth = 0.5,
                               title_suffix = " (single pixel)")
ggsave(file.path(output_dir, "NDVI_spiral_single_pixel.png"),
       p_single, width = 6, height = 6, dpi = 300, bg = "white")
message("\u2705 Single-pixel spiral saved")

# ── All pixels ───────────────────────────────────────────
all_paths <- interpolate_arcs(df_long, n_interp = 5)
p_all <- build_spiral_plot(all_paths, grid, yr_range)
ggsave(file.path(output_dir, "NDVI_spiral_Wild_Ennerdale.png"),
       p_all, width = 6, height = 6, dpi = 300, bg = "white")
message("\u2705 Full spiral saved")
