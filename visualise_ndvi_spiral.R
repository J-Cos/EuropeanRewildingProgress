# ==============================================================================
# NDVI Seasonality Spiral — Wild Ennerdale (demo site)
# Polar plot: angle = day-of-year, radius = NDVI, colour = year
# Each pixel is a semi-transparent line; colour gradient from first to last year
# PNAS style
# ==============================================================================

library(terra)
library(tidyverse)
library(cowplot)

# ── Config ───────────────────────────────────────────────
input_file <- "Outputs/GEE_MK_250m/NDVI_TimeSeries_250m_Wild_Ennerdale.tif"
output_dir <- "Outputs/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load raster and parse dates from band names ──────────
r <- rast(input_file)
NAflag(r) <- -9999

band_names <- names(r)
# Band names: NDVI_YYYYMMDD → parse dates
date_strings <- gsub("^NDVI_", "", band_names)
dates <- as.Date(date_strings, format = "%Y%m%d")

# ── Extract pixel values → long format ───────────────────
vals <- values(r)  # matrix: rows = pixels, cols = dates
colnames(vals) <- as.character(dates)

# Remove fully-NA pixels
valid_px <- rowSums(!is.na(vals)) > 0
vals <- vals[valid_px, , drop = FALSE]

# Add pixel IDs and pivot long
df <- as.data.frame(vals)
df$pixel_id <- seq_len(nrow(df))

df_long <- df %>%
  pivot_longer(
    cols = -pixel_id,
    names_to = "date_str",
    values_to = "ndvi"
  ) %>%
  mutate(
    date = as.Date(date_str),
    year = lubridate::year(date),
    doy  = lubridate::yday(date)
  ) %>%
  filter(!is.na(ndvi), ndvi >= 0, ndvi <= 1) %>%
  select(pixel_id, date, year, doy, ndvi)

# ── Month labels for polar axis ──────────────────────────
month_breaks <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
month_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# ── Year range for colour scale ──────────────────────────
yr_range <- range(df_long$year)

# ── Spiral plot ──────────────────────────────────────────
p <- ggplot(df_long, aes(
    x = doy,
    y = ndvi,
    group = interaction(pixel_id, year),
    colour = year
  )) +
  geom_line(alpha = 0.015, linewidth = 0.25) +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_x_continuous(
    breaks = month_breaks,
    labels = month_labels,
    limits = c(1, 366),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0.2, 0.8, by = 0.2),
    expand = c(0, 0)
  ) +
  scale_colour_viridis_c(
    name = NULL,
    option = "inferno",
    begin = 0.15,
    end = 0.9,
    breaks = pretty(yr_range[1]:yr_range[2], n = 5),
    guide = guide_colourbar(
      barwidth = unit(0.3, "cm"),
      barheight = unit(4, "cm"),
      ticks = FALSE,
      title.position = "top"
    )
  ) +
  labs(
    title = "Wild Ennerdale — NDVI Seasonality",
    subtitle = paste0(
      "Each line = one pixel\u2013year (",
      yr_range[1], "\u2013", yr_range[2], ")"
    )
  ) +
  theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(
    # PNAS-style clean layout
    plot.title       = element_text(size = 10, face = "bold", hjust = 0.5,
                                     margin = margin(b = 2)),
    plot.subtitle    = element_text(size = 7.5, hjust = 0.5, colour = "grey40",
                                     margin = margin(b = 8)),
    axis.text.x      = element_text(size = 7, colour = "grey30"),
    axis.text.y      = element_text(size = 6, colour = "grey50"),
    axis.title       = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position  = c(0.92, 0.15),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 7),
    plot.background  = element_rect(fill = "white", colour = NA),
    plot.margin      = margin(10, 10, 10, 10)
  )

# ── Save ─────────────────────────────────────────────────
out_file <- file.path(output_dir, "NDVI_spiral_Wild_Ennerdale.png")
ggsave(out_file, p, width = 6, height = 6, dpi = 300, bg = "white")
message(paste("\u2705 Spiral plot saved to", out_file))
