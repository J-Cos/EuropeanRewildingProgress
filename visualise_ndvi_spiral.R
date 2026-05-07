# ==============================================================================
# NDVI Seasonality Spiral — Wild Ennerdale (demo site)
# Polar plot: angle = day-of-year, radius = NDVI, colour = year
# Each pixel is a single continuous path; segments wrap at the Jan boundary
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
date_strings <- gsub("^NDVI_", "", band_names)
dates <- as.Date(date_strings, format = "%Y%m%d")

# ── Extract pixel values → long format ───────────────────
vals <- values(r)
colnames(vals) <- as.character(dates)

# Remove fully-NA pixels
valid_px <- rowSums(!is.na(vals)) > 0
vals <- vals[valid_px, , drop = FALSE]

df <- as.data.frame(vals)
df$pixel_id <- seq_len(nrow(df))

df_long <- df %>%
  pivot_longer(cols = -pixel_id, names_to = "date_str", values_to = "ndvi") %>%
  mutate(
    date = as.Date(date_str),
    year = lubridate::year(date),
    doy  = lubridate::yday(date)
  ) %>%
  filter(!is.na(ndvi), ndvi >= 0, ndvi <= 1) %>%
  select(pixel_id, date, year, doy, ndvi) %>%
  arrange(pixel_id, date)

# ── Build segments between consecutive observations ──────
# Each pixel is one continuous time series; we draw segment-by-segment
# so that Dec→Jan transitions wrap correctly around the polar axis.
df_segs <- df_long %>%
  group_by(pixel_id) %>%
  mutate(
    doy_end  = lead(doy),
    ndvi_end = lead(ndvi),
    year_end = lead(year)
  ) %>%
  filter(!is.na(doy_end)) %>%
  ungroup()

# Separate normal segments from those that cross the year boundary
normal <- df_segs %>% filter(doy_end >= doy)
wrapping <- df_segs %>% filter(doy_end < doy)

# Split wrapping segments at the boundary (doy = 366 / 0)
# Interpolate NDVI at the wrap point
wrap_to_dec <- wrapping %>%
  mutate(
    # fractional distance from start to the 366 boundary
    frac = (366 - doy) / (366 - doy + doy_end),
    ndvi_boundary = ndvi + frac * (ndvi_end - ndvi),
    doy_end = 366,
    ndvi_end = ndvi_boundary
  ) %>%
  select(pixel_id, year, doy, ndvi, doy_end, ndvi_end)

wrap_from_jan <- wrapping %>%
  mutate(
    frac = (366 - doy) / (366 - doy + doy_end),
    ndvi_boundary = ndvi + frac * (ndvi_end - ndvi),
    doy = 0,
    ndvi = ndvi_boundary,
    year = year_end
  ) %>%
  select(pixel_id, year, doy, ndvi, doy_end, ndvi_end)

all_segs <- bind_rows(
  normal %>% select(pixel_id, year, doy, ndvi, doy_end, ndvi_end),
  wrap_to_dec,
  wrap_from_jan
)

# ── Plot config ──────────────────────────────────────────
month_breaks <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
month_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
yr_range <- range(df_long$year)

# ── Spiral plot ──────────────────────────────────────────
p <- ggplot(all_segs, aes(
    x = doy, xend = doy_end,
    y = ndvi, yend = ndvi_end,
    colour = year
  )) +
  geom_segment(alpha = 0.012, linewidth = 0.2) +
  coord_polar(start = -pi / 2, direction = 1) +
  scale_x_continuous(
    breaks = month_breaks,
    labels = month_labels,
    limits = c(0, 366),
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
    title = "Wild Ennerdale \u2014 NDVI Seasonality",
    subtitle = paste0(
      "Each line = one pixel (",
      yr_range[1], "\u2013", yr_range[2], ")"
    )
  ) +
  theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(
    plot.title       = element_text(size = 10, face = "bold", hjust = 0.5,
                                     margin = margin(b = 2)),
    plot.subtitle    = element_text(size = 7.5, hjust = 0.5, colour = "grey40",
                                     margin = margin(b = 8)),
    axis.text.x      = element_text(size = 7, colour = "grey30"),
    axis.text.y      = element_text(size = 6, colour = "grey50"),
    axis.title       = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position.inside = c(0.92, 0.15),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 7),
    plot.background  = element_rect(fill = "white", colour = NA),
    plot.margin      = margin(10, 10, 10, 10)
  )

# ── Save ─────────────────────────────────────────────────
out_file <- file.path(output_dir, "NDVI_spiral_Wild_Ennerdale.png")
ggsave(out_file, p, width = 6, height = 6, dpi = 300, bg = "white")
message(paste("\u2705 Spiral plot saved to", out_file))
