#!/usr/bin/env Rscript

# the canonical 217-economy sample. source this
#
# 8 world bank aggregates (S1-S4, V1-V4) pool economies already counted here, 
# so they  double-count and are dropped. the 24 non-sovereign units are KEPT:
# dropping them is not neutral, it raises developed adherence from 41.2 to 42.9 and weakens 
# the paper's third claim

PANEL_AGGREGATES <- c("S1", "S2", "S3", "S4", "V1", "V2", "V3", "V4")

load_panel <- function(drop_aggregates = TRUE) {
  base <- here::here()
  d <- utils::read.csv(file.path(base, "data/magic_squares_data_clean.csv"))
  # Namibia iso-2 is the string "NA"; read.csv makes it missing. repair here, not
  # via na.strings, which would parse every numeric column as text
  d$iso2c[d$country == "Namibia"] <- "NA"
  if (drop_aggregates) d <- d[!d$iso2c %in% PANEL_AGGREGATES, ]
  d <- d[order(d$country, d$year), ]
  rownames(d) <- NULL
  d
}

panel_note <- function(d)
  sprintf("%d economies, %d country-years, %d-%d",
          length(unique(d$country)), nrow(d), min(d$year), max(d$year))
