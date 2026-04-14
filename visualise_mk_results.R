# ==============================================================================
# Multi-panel figure: Sen's Slope (INDVI) — tidyterra + cowplot
# Equal square panels, exact 1 km scale bars, shared legend
# ==============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)

input_dir  <- "Outputs/GEE_MK_250m"
output_dir <- "Outputs/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Site metadata ────────────────────────────────────────
meta <- data.frame(
  site_id = c(
    "Alladale", "Border_Meuse", "Central_Apennines", "Chernobyl_Exclusion_Zone",
    "Creag_Meagaidh", "Doberitzer_Heide", "Greater_Coa_Valley", "Naliboki_Forest",
    "Oostvaardersplassen", "Pentezug", "Swiss_National_Park",
    "Knepp", "Wild_Ennerdale", "Tarutino_Steppe",
    "Dundreggan", "Kempen_Broek",
    "Southern_Carpathians", "Velebit_Mountains", "Vanatori_Neamt_Nature_Park",
    "Danube_Delta"
  ),
  start = c(
    2000, 2000, 2000, 2000,
    2000, 2000, 2000, 2000,
    2000, 2000, 2000,
    2001, 2003, 2004,
    2008, 2010,
    2011, 2012, 2012,
    2013
  ),
  end = 2025,
  stringsAsFactors = FALSE
)
meta <- meta[order(meta$start, meta$site_id), ]

# ── Load layers ──────────────────────────────────────────
slopes <- list()
titles <- character()
for (i in seq_len(nrow(meta))) {
  sid <- meta$site_id[i]
  fname <- file.path(input_dir, paste0("MK_Stacked_250m_", sid, ".tif"))
  if (!file.exists(fname)) { message("  Skipping: ", sid); next }
  r <- rast(fname)
  NAflag(r) <- -9999
  slopes[[length(slopes) + 1]] <- r[["SenSlope_INDVI"]]
  label <- gsub("_", " ", sid)
  titles <- c(titles, paste0(label, " (", meta$start[i], "\u2013", meta$end[i], ")"))
}

# ── Shared colour limits ─────────────────────────────────
all_vals <- unlist(lapply(slopes, function(r) values(r, na.rm = TRUE)))
zlim <- c(-1, 1) * max(abs(quantile(all_vals, c(0.02, 0.98))), na.rm = TRUE)

shared_scale <- scale_fill_gradientn(
  colours = hcl.colors(50, "RdYlGn"),
  limits = zlim,
  na.value = "transparent",
  name = expression("Sen's Slope (INDVI yr"^{-1}*")")
)

# ── Build individual panels (square extent, 1km bar) ─────
make_panel <- function(r, title) {
  e <- ext(r)
  lat_mid <- (e$ymin + e$ymax) / 2
  cos_lat <- cos(lat_mid * pi / 180)

  # Convert extent to approximate metres, find the longer side, pad to square
  w_m <- (e$xmax - e$xmin) * 111320 * cos_lat
  h_m <- (e$ymax - e$ymin) * 110540
  side_m <- max(w_m, h_m) * 1.15  # 15% padding

  cx <- (e$xmin + e$xmax) / 2
  cy <- (e$ymin + e$ymax) / 2
  half_x_deg <- (side_m / 2) / (111320 * cos_lat)
  half_y_deg <- (side_m / 2) / 110540

  xlims <- c(cx - half_x_deg, cx + half_x_deg)
  ylims <- c(cy - half_y_deg, cy + half_y_deg)

  # 1 km bar in degrees
  bar_deg <- 1000 / (111320 * cos_lat)

  # Scale bar position: bottom-left of the padded square
  sb_x <- xlims[1] + (xlims[2] - xlims[1]) * 0.05
  sb_y <- ylims[1] + (ylims[2] - ylims[1]) * 0.05

  ggplot() +
    geom_spatraster(data = r) +
    shared_scale +
    annotate("segment",
             x = sb_x, xend = sb_x + bar_deg,
             y = sb_y, yend = sb_y,
             linewidth = 0.8, colour = "black") +
    coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 6.5, margin = margin(b = 2)),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2)
    )
}

panels <- mapply(make_panel, slopes, titles, SIMPLIFY = FALSE)

# Pad to 20 with blank plots
while (length(panels) < 20) {
  panels[[length(panels) + 1]] <- ggplot() + theme_void()
}

# ── Extract shared legend ────────────────────────────────
legend_src <- ggplot() +
  geom_spatraster(data = slopes[[1]]) +
  shared_scale +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key.height = unit(2.5, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  )
shared_legend <- get_legend(legend_src)

# ── Scale bar annotation for legend column ───────────────
bar_label <- ggdraw() +
  draw_line(x = c(0.2, 0.7), y = c(0.6, 0.6), size = 1, color = "black") +
  draw_label("= 1 km", x = 0.5, y = 0.35, size = 7, fontface = "plain")

legend_col <- plot_grid(shared_legend, bar_label, ncol = 1, rel_heights = c(1, 0.08))

# ── Assemble ─────────────────────────────────────────────
grid <- plot_grid(plotlist = panels, ncol = 5)
final <- plot_grid(grid, legend_col, ncol = 2, rel_widths = c(1, 0.08))

out_file <- file.path(output_dir, "SenSlope_INDVI_multipanel.png")
ggsave(out_file, final, width = 15, height = 12, dpi = 200, bg = "white")
message(paste("\u2705 Multi-panel figure saved to", out_file))
