library(testthat)
library(dplyr)
library(terra)
library(lubridate)

source("visualise_ndvi_spiral.R")

test_that("polar_to_cartesian computes correctly without jitter", {
  df <- tibble(pixel_id = 1, doy = 365.25/4, ndvi = 1)
  res <- polar_to_cartesian(df)
  
  # doy = 91.3125 (quarter circle), angle should be -pi/2 - pi/2 = -pi
  expect_equal(res$cx, -1, tolerance = 1e-5)
  expect_equal(res$cy, 0, tolerance = 1e-5)
})

test_that("polar_to_cartesian with jitter changes doy slightly", {
  set.seed(42)
  df <- tibble(pixel_id = 1:100, doy = rep(100, 100), ndvi = rep(0.5, 100))
  res <- polar_to_cartesian(df, jitter_days = 8)
  
  expect_true(any(res$doy_j != 100))
  expect_true(all(abs(res$doy_j - 100) <= 8))
})

test_that("interpolate_arcs breaks on large gaps and interpolates properly", {
  df <- tibble(
    pixel_id = c(1, 1, 1),
    date = as.Date(c("2010-01-01", "2010-01-17", "2010-06-01")),
    year = c(2010, 2010, 2010),
    doy = yday(date),
    ndvi = c(0.2, 0.4, 0.8)
  )
  
  res <- interpolate_arcs(df, n_interp = 3, max_gap_days = 50)
  
  # The gap between 2010-01-17 and 2010-06-01 is > 50 days (135 days)
  # It should result in an NA break
  na_rows <- res %>% filter(is.na(ndvi))
  expect_equal(nrow(na_rows), 1)
  
  # + 2 endpoints = 5 points. The NA break also emits the segment start point
  # again before the NA, so there are 6 valid rows.
  valid_rows <- res %>% filter(!is.na(ndvi))
  expect_equal(nrow(valid_rows), 6)
})

test_that("make_grid returns expected structures", {
  g <- make_grid()
  expect_named(g, c("rings", "spokes", "ring_labels"))
  expect_equal(nrow(g$spokes), 12) # 12 months
})
