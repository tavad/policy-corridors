#!/usr/bin/env Rscript
# figures 1-4. every value is read from the table files at run time, so a figure
# cannot drift from the table it illustrates.
# matplotlib blue/orange: dE 24.6 under protanopia, and they differ in lightness so
# they survive a greyscale print; shape and fill repeat the split, colour is never alone

options(stringsAsFactors = FALSE)
BASE <- here::here()
OUT  <- file.path(BASE, "plots")
dir.create(OUT, showWarnings = FALSE)

BLUE <- "#1f77b4"; ORANGE <- "#ff7f0e"
INK  <- "#1a1a1a"; MUTED <- "#5a5a5a"; GRID <- "#d8d8d8"; PALE <- "#c3ced6"

HAVE_RAGG <- requireNamespace("ragg", quietly = TRUE)
open_png <- function(name, w, h) {
  path <- file.path(OUT, name)
  if (HAVE_RAGG) ragg::agg_png(path, width = w, height = h, units = "in", res = 300,
                               background = "white", pointsize = 9)
  else png(path, width = w, height = h, units = "in", res = 300, type = "cairo",
           bg = "white", pointsize = 9)
  path
}

LH_PT <- 10.5; PAD_PT <- 9
cap <- function(txt) {
  ln <- strsplit(txt, "\n", fixed = TRUE)[[1]]; n <- length(ln)
  hpt <- par("din")[2] * 72
  x <- grconvertX(0.012, from = "ndc", to = "user")
  for (i in seq_along(ln))
    text(x, grconvertY((PAD_PT + (n - i) * LH_PT) / hpt, from = "ndc", to = "user"),
         ln[i], adj = c(0, 0), cex = 0.78, col = MUTED, xpd = NA)
}
title_left <- function(m, cex = 1.17, y = 0.965)
  text(grconvertX(0.012, "ndc", "user"), grconvertY(y, "ndc", "user"),
       m, adj = c(0, 1), cex = cex, col = INK, xpd = NA)
subtitle_left <- function(s, cex = 0.85, y = 0.905)
  text(grconvertX(0.012, "ndc", "user"), grconvertY(y, "ndc", "user"),
       s, adj = c(0, 1), cex = cex, col = MUTED, xpd = NA)
overlay <- function() { par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
                        plot.new(); plot.window(c(0, 1), c(0, 1)) }

IND_LAB <- c(real_interest = "Real interest rate", inflation = "Inflation",
             broad_money_gdp = "Broad money / GDP", trade_balance = "Trade balance",
             unemployment = "Unemployment")
ORDER <- names(IND_LAB)

# F1: adherence by indicator and development group, 2008-2023
# same file Table 3 prints from, so figure and table cannot disagree; the earlier
# version read the 225-entity result and put four bars up to 2.7 points off
t3  <- read.csv(file.path(BASE, "outputs/T3_adherence.csv"))
dev <- setNames(t3$dev_pct, t3$indicator)[ORDER]
tra <- setNames(t3$tra_pct, t3$indicator)[ORDER]

open_png("F1_adherence.png", 7.2, 5.3)
par(mar = c(6.4, 12.0, 4.4, 1.4), family = "sans", las = 1)
plot.new(); plot.window(xlim = c(0, 70), ylim = c(0.4, length(ORDER) + 0.6))
abline(v = seq(0, 70, 10), col = GRID, lwd = 0.7)
abline(v = 50, col = MUTED, lwd = 1.1, lty = 2)
h <- 0.30
for (i in seq_along(ORDER)) {
  y <- length(ORDER) - i + 1
  rect(0, y + 0.04, dev[i], y + 0.04 + h, col = BLUE,   border = NA)
  rect(0, y - 0.04 - h, tra[i], y - 0.04, col = ORANGE, border = NA,
       density = 18, angle = 45)
  rect(0, y - 0.04 - h, tra[i], y - 0.04, col = NA, border = ORANGE, lwd = 1.2)
  text(dev[i] + 1.2, y + 0.04 + h / 2, sprintf("%.1f", dev[i]), adj = c(0, 0.5), cex = 0.82, col = INK)
  text(tra[i] + 1.2, y - 0.04 - h / 2, sprintf("%.1f", tra[i]), adj = c(0, 0.5), cex = 0.82, col = INK)
  mtext(IND_LAB[ORDER[i]], side = 2, at = y, line = 0.6, adj = 1, cex = 0.95, col = INK)
}
axis(1, at = seq(0, 70, 10), col = MUTED, col.axis = MUTED, cex.axis = 0.85, lwd = 0.8)
mtext("Share of country-years inside the corridor (%)", side = 1, line = 2.4, cex = 0.95, col = INK)
text(50, length(ORDER) + 0.55, "50% (majority)", adj = c(0.5, 0), cex = 0.78, col = MUTED, xpd = NA)
overlay()
title_left("Adherence is a minority condition")
subtitle_left("Share of country-years inside each corridor, 2008-2023", y = 0.898)
legend(0.60, 0.99, c("Developed", "Transforming"), fill = c(BLUE, ORANGE),
       density = c(NA, 18), angle = c(NA, 45), border = c(NA, ORANGE),
       bty = "n", cex = 0.85, text.col = INK, horiz = TRUE, xpd = NA)
cap(paste("Source: outputs/T3_adherence.csv, the same computation Table 3 is printed from: 217 economies, 2008-2023,",
          "country-years where the indicator and GDP growth are both observed. Only the developed trade balance and",
          "unemployment corridors are met by a majority.",
          sep = "\n"))
invisible(dev.off())

# F2: adherence 1960-2023, all countries against a balanced set
ay <- read.csv(file.path(BASE, "data/historical_20260816/adherence_by_year.csv"))

open_png("F2_adherence_1960_2023.png", 7.2, 6.0)
layout(matrix(1:2, ncol = 1), heights = c(2.5, 1))
par(oma = c(5.2, 0, 0, 0), mar = c(0.4, 4.6, 4.2, 1.4), family = "sans", las = 1)
plot.new(); plot.window(xlim = range(ay$year), ylim = c(0.5, 1.75))
abline(h = seq(0.6, 1.7, 0.2), col = GRID, lwd = 0.7)
lines(ay$year, ay$mean_K3_all, col = PALE, lwd = 2.4)
lines(ay$year, ay$mean_K3_balanced, col = BLUE, lwd = 2.4)
axis(2, at = seq(0.6, 1.6, 0.2), col = MUTED, col.axis = MUTED, cex.axis = 0.85, lwd = 0.8)
mtext("Corridors met, mean of 3", side = 2, line = 3.0, las = 0, cex = 0.95, col = INK)
# the three largest annual falls, marked because Table 6 attributes each to a corridor
for (yy in c(1974, 2008, 2022)) {
  g <- ay[ay$year == yy, ]
  points(g$year, g$mean_K3_all, pch = 21, col = MUTED, bg = "white", cex = 1.1)
  text(g$year, g$mean_K3_all - 0.055, sprintf("%d: %.2f", yy, g$mean_K3_all),
       adj = c(0.5, 1), cex = 0.75, col = MUTED)
}
legend("topleft", c("54 economies scored throughout", "all economies scored that year"),
       col = c(BLUE, PALE), lwd = 2.4, bty = "n", cex = 0.85, text.col = INK)
par(mar = c(2.9, 4.6, 0.6, 1.4))
plot.new(); plot.window(xlim = range(ay$year), ylim = c(0, max(ay$n_countries_all) * 1.08))
abline(h = seq(0, 140, 40), col = GRID, lwd = 0.7)
polygon(c(ay$year, rev(ay$year)), c(ay$n_countries_all, rep(0, nrow(ay))),
        col = "#e8edf1", border = NA)
lines(ay$year, ay$n_countries_all, col = MUTED, lwd = 1.4)
axis(1, at = seq(1960, 2020, 10), col = MUTED, col.axis = MUTED, cex.axis = 0.85, lwd = 0.8)
axis(2, at = seq(0, 120, 40), col = MUTED, col.axis = MUTED, cex.axis = 0.85, lwd = 0.8)
mtext("Countries scored", side = 2, line = 3.0, las = 0, cex = 0.95, col = INK)
overlay()
title_left("Adherence fell through the 1970s, then rose after 1980")
subtitle_left("Mean number of the three corridors met. The rise dates from 1980; the labelled years are the three largest falls.")
cap(paste("Source: data/historical_20260816/adherence_by_year.csv. Inflation, broad money/GDP and the trade balance, the three series",
          "reaching back to 1960. Fitted over the whole period the full sample is flat (+0.018 corridors per decade, CI -0.006 to",
          "+0.042); fitted from 1980 it rises (+0.096, +0.068 to +0.123). The 54 economies scored throughout rise on both windows,",
          "+0.064 and +0.147. Later entrants are less compliant, so the growing sample understates the rise within it.",
          sep = "\n"))
invisible(dev.off())

# F3: which bounds behave like corridors
# volatility, not the success rate: on success the developed trade balance looks
# monotone (39.1/43.6/47.9) but on volatility it is the corridor shape (2.92/2.46/4.29)
asym <- read.csv(file.path(BASE, "data/necessity_20260816/asymmetry_by_side.csv"))
SIDES <- c("below", "inside", "above")
GRPS  <- c("Developed", "Developing")
GLAB  <- c(Developed = "Developed", Developing = "Transforming")

open_png("F3_corridor_test.png", 7.8, 5.4)
layout(matrix(1:10, nrow = 2, byrow = TRUE))
par(oma = c(6.2, 4.4, 5.0, 0.8), mar = c(0.7, 0.6, 1.8, 0.4), family = "sans", las = 1)
for (g in GRPS) for (ind in ORDER) {
  z <- asym[asym$group == g & asym$indicator == ind, ]
  v <- setNames(rep(NA_real_, 3), SIDES)
  if (nrow(z)) v[z$side] <- z$sd_growth_5y
  plot.new(); plot.window(xlim = c(0.6, 3.4), ylim = c(2.0, 6.6))
  abline(h = 2:6, col = GRID, lwd = 0.7)
  lines(1:3, v, col = PALE, lwd = 1.8)
  best <- which.min(v)
  for (k in 1:3) {
    ins <- SIDES[k] == "inside"
    points(k, v[k], pch = if (ins) 19 else 17,
           col = if (ins) BLUE else ORANGE, cex = if (ins) 1.5 else 1.15)
    points(k, v[k], pch = if (ins) 21 else 24, col = "white", bg = NA,
           cex = if (ins) 1.5 else 1.15, lwd = 0.9)
  }
  if (length(best) && SIDES[best] == "inside")
    text(2, 2.12, "steadiest inside", cex = 0.68, col = BLUE, adj = c(0.5, 0.5))
  if (g == "Developed") mtext(IND_LAB[ind], side = 3, line = 0.35, cex = 0.76, col = INK)
  if (ind == ORDER[1]) {
    axis(2, at = 2:6, col = MUTED, col.axis = MUTED, cex.axis = 0.78, lwd = 0.8)
    text(0.66, 6.42, GLAB[g], adj = c(0, 0.5), cex = 0.95, col = INK, font = 2)
  }
  if (g == "Developing") axis(1, at = 1:3, labels = SIDES, col = MUTED, col.axis = MUTED,
                              cex.axis = 0.76, lwd = 0.8)
}
overlay()
title_left("Only some bounds behave like corridors")
subtitle_left("Volatility of GDP growth over five years, by side of each bound. A corridor is lowest in the middle.", y = 0.925)
mtext("Standard deviation of growth (pp)", side = 2, line = 2.6, outer = TRUE, las = 0,
      cex = 0.9, col = INK)
cap(paste("Source: data/necessity_20260816/asymmetry_by_side.csv. Three of the five developed bounds are steadiest inside and two of the",
          "transforming ones are. Where they fail it is the floor: developed inflation and broad money are steadiest below their lower bound,",
          "and transforming unemployment is monotone, so its 5% floor has no support. Success rates and worst-year depth are in Table 5.",
          sep = "\n"))
invisible(dev.off())

# F4: Korea and China, compliance and growth move in opposite directions
# the 164-take-off event study is a sentence in the text, not drawn: its bootstrap band
# overlaps the panel level, so a figure would overstate it; the per-economy path is legible
cdec <- read.csv(file.path(BASE, "data/historical_20260816/case_decade.csv"))
et   <- read.csv(file.path(BASE, "data/historical_20260816/event_time_series.csv"))
PANEL_REF <- et$panel_ref[1]

CASE <- c("Korea, Rep.", "China"); CLAB <- c("Korea", "China")
CCOL <- c(BLUE, ORANGE); CPCH <- c(19, 17)

open_png("F4_korea_china.png", 7.4, 5.4)
par(oma = c(5.0, 0, 0, 0), mar = c(4.2, 4.8, 4.6, 1.6), family = "sans", las = 1)
plot.new(); plot.window(xlim = c(0, 2.6), ylim = c(2.5, 11.4))
abline(h = seq(3, 11, 2), col = GRID, lwd = 0.7)
abline(v = seq(0, 2.5, 0.5), col = GRID, lwd = 0.7)
abline(v = PANEL_REF, col = MUTED, lwd = 1.0, lty = 2)
text(PANEL_REF - 0.05, 3.9, sprintf("panel-wide mean, %.2f", PANEL_REF),
     adj = c(0, 0.5), cex = 0.76, col = MUTED, srt = 90)

for (i in seq_along(CASE)) {
  z <- cdec[cdec$country == CASE[i] & cdec$n_years >= 5, ]
  z <- z[order(z$decade), ]
  np <- nrow(z)
  lines(z$K3, z$gdp_growth, col = CCOL[i], lwd = 2.0)
  arrows(z$K3[np - 1], z$gdp_growth[np - 1], z$K3[np], z$gdp_growth[np],
         col = CCOL[i], lwd = 2.0, length = 0.11)
  points(z$K3, z$gdp_growth, pch = CPCH[i], col = CCOL[i], cex = 1.35)
  points(z$K3, z$gdp_growth, pch = if (CPCH[i] == 19) 21 else 24, col = "white",
         bg = NA, cex = 1.35, lwd = 1.0)
  # decade labels sit above or below the marker depending on where the point is
  # in the country's own range, so a label never lands on the line
  off <- ifelse(z$gdp_growth >= mean(range(z$gdp_growth)), 0.36, -0.36)
  text(z$K3, z$gdp_growth + off, paste0(substr(z$decade, 3, 4), "s"),
       cex = 0.78, col = CCOL[i], adj = c(0.5, 0.5))
}

axis(1, at = seq(0, 2.5, 0.5), col = MUTED, col.axis = MUTED, cex.axis = 0.85, lwd = 0.8)
axis(2, at = seq(3, 11, 2), col = MUTED, col.axis = MUTED, cex.axis = 0.85, lwd = 0.8)
mtext("Corridors met, decade mean (of 3)", side = 1, line = 2.6, cex = 0.95, col = INK)
mtext("GDP growth, decade mean (%)", side = 2, line = 3.2, las = 0, cex = 0.95, col = INK)
legend("topright", CLAB, col = CCOL, pch = CPCH, lwd = 2.0, bty = "n",
       cex = 0.88, text.col = INK)
overlay()
title_left("Compliance rose as growth fell")
subtitle_left("Corridor compliance against GDP growth, decade by decade. Arrows run forward in time.")
cap(paste("Source: data/historical_20260816/case_decade.csv, built from data/wdi_long_1960_2023.csv. Korea moves right and down: its",
          "lowest compliance came in the 9-10% decades and its highest in the 2000s, when growth had fallen to 4.9%. China's",
          "fastest decades, the 1990s and 2000s, are its least compliant. Partial decades of fewer than five scored years are",
          "omitted. Across all 164 take-offs identified panel-wide, compliance in the run-up is 0.95 against 1.08 elsewhere.",
          sep = "\n"))
invisible(dev.off())

cat("written:\n"); for (f in list.files(OUT, pattern = "^F[0-9].*png$")) cat("  plots/", f, "\n", sep = "")
