#!/usr/bin/env Rscript

# does the corridor framework describe an earlier monetary order? long panel, 1960-2023.
# take-offs are identified by a rule on every country, not on the famous cases:
# selecting on the outcome is the error the manuscript concedes for china/singapore/korea.
# unemployment is absent in wdi before 1991 and the real rate is thin before 1976, so the
# historical score uses the three back-reaching series: inflation, broad money/gdp, trade balance

options(width = 175, stringsAsFactors = FALSE)
suppressWarnings(RNGversion(getRversion()))

BASE <- here::here()
OUT  <- file.path(BASE, "data/historical_20260816")
dir.create(OUT, showWarnings = FALSE)
SEED <- 42L; B_BOOT <- 20000L
set.seed(SEED)

logcon <- file(file.path(OUT, "historical_run_log.txt"), open = "wt")
sink(logcon, split = TRUE); on.exit({ sink(); close(logcon) }, add = TRUE)
msg <- function(...) cat(..., "\n", sep = "")
hdr <- function(s) cat("\n", strrep("=", 104), "\n", s, "\n", strrep("=", 104), "\n", sep = "")

hdr("historical_test.R -- corridor compliance and growth take-offs, 1960-2023")
msg("run: ", format(Sys.time()), "  seed=", SEED)

d <- read.csv(file.path(BASE, "data/wdi_long_1960_2023.csv"))
d <- d[order(d$country, d$year), ]
msg("rows=", nrow(d), "  countries=", length(unique(d$country)),
    "  years=", min(d$year), "-", max(d$year))

# status assigned per country-year from gdp per capita (12,696 high-income line), so a
# country changes stage over the period instead of being averaged into one
HI_LINE <- 12696
d$status <- ifelse(is.na(d$gdp_pc), NA,
                   ifelse(d$gdp_pc >= HI_LINE, "Developed", "Developing"))

CORR <- list(
  Developed  = list(inflation = c(1, 3), broad_money_gdp = c(70, 120), trade_balance = c(-3, 6)),
  Developing = list(inflation = c(2, 5), broad_money_gdp = c(40, 80),  trade_balance = c(-8, 3)))
SET3 <- c("inflation", "broad_money_gdp", "trade_balance")

for (ind in SET3) {
  lo <- sapply(d$status, function(g) if (is.na(g)) NA else CORR[[g]][[ind]][1])
  hi <- sapply(d$status, function(g) if (is.na(g)) NA else CORR[[g]][[ind]][2])
  d[[paste0("in_", ind)]] <- ifelse(is.na(d[[ind]]) | is.na(lo), NA,
                                    as.numeric(d[[ind]] >= lo & d[[ind]] < hi))
}
d$K3 <- ifelse(rowSums(is.na(d[, paste0("in_", SET3)])) > 0, NA,
               rowSums(d[, paste0("in_", SET3)]))
msg("scored (K3) country-years: ", sum(!is.na(d$K3)),
    "  countries: ", length(unique(d$country[!is.na(d$K3)])))

key <- paste(d$country, d$year)
lagv  <- function(v, h) v[match(paste(d$country, d$year - h), key)]
leadv <- function(v, h) v[match(paste(d$country, d$year + h), key)]

hdr("BLOCK 1 -- adherence and its association with growth, by era")

d$era <- cut(d$year, breaks = c(1959, 1969, 1979, 1989, 1999, 2007, 2023),
             labels = c("1960s", "1970s", "1980s", "1990s", "2000-07", "2008-23"))

FW <- sapply(1:5, function(h) leadv(d$gdp_growth, h))
ok5 <- rowSums(is.na(FW)) == 0
d$g_fwd5   <- ifelse(ok5, rowMeans(FW), NA)
d$sd_fwd5  <- ifelse(ok5, apply(FW, 1, sd), NA)

# country-cluster bootstrap on the high-minus-low gap, so a thin era cannot look confident
gap_ci <- function(y, hi, lo, cid) {
  ok <- !is.na(y) & (hi | lo); y <- y[ok]; hi <- hi[ok]; lo <- lo[ok]; cid <- cid[ok]
  if (sum(hi) < 5 || sum(lo) < 5) return(c(NA, NA))
  cf <- factor(cid); i <- as.integer(cf); C <- nlevels(cf)
  Sh <- as.vector(rowsum(y * hi, i, reorder = TRUE)); Nh <- as.vector(rowsum(as.numeric(hi), i, reorder = TRUE))
  Sl <- as.vector(rowsum(y * lo, i, reorder = TRUE)); Nl <- as.vector(rowsum(as.numeric(lo), i, reorder = TRUE))
  M <- stats::rmultinom(B_BOOT, size = C, prob = rep(1 / C, C))
  bt <- ifelse(as.vector(crossprod(M, Nh)) > 0 & as.vector(crossprod(M, Nl)) > 0,
               as.vector(crossprod(M, Sh)) / as.vector(crossprod(M, Nh)) -
               as.vector(crossprod(M, Sl)) / as.vector(crossprod(M, Nl)), NA_real_)
  unname(quantile(bt, c(.025, .975), na.rm = TRUE))
}

rows <- list()
for (er in levels(d$era)) for (st in c("Developed", "Developing")) {
  s <- which(d$era == er & d$status == st & !is.na(d$K3))
  if (length(s) < 40) next
  hiK <- d$K3[s] >= 2; loK <- d$K3[s] <= 1
  gci <- gap_ci(d$g_fwd5[s],  hiK, loK, d$country[s])
  sci <- gap_ci(d$sd_fwd5[s], hiK, loK, d$country[s])
  rows[[length(rows) + 1]] <- data.frame(
    era = er, status = st, n = length(s), n_countries = length(unique(d$country[s])),
    pct_adherent_all3 = 100 * mean(d$K3[s] == 3),
    mean_K = mean(d$K3[s]),
    growth_fwd5_hiK = mean(d$g_fwd5[s][hiK], na.rm = TRUE),
    growth_fwd5_loK = mean(d$g_fwd5[s][loK], na.rm = TRUE),
    sd_fwd5_hiK = mean(d$sd_fwd5[s][hiK], na.rm = TRUE),
    sd_fwd5_loK = mean(d$sd_fwd5[s][loK], na.rm = TRUE),
    growth_gap_lo = gci[1], growth_gap_hi = gci[2],
    sd_gap_lo = sci[1], sd_gap_hi = sci[2])
}
ERA <- do.call(rbind, rows)
ERA$growth_gap <- ERA$growth_fwd5_hiK - ERA$growth_fwd5_loK
ERA$sd_gap     <- ERA$sd_fwd5_hiK - ERA$sd_fwd5_loK
print(ERA[, c("era", "status", "n", "n_countries", "pct_adherent_all3", "mean_K",
              "growth_gap", "growth_gap_lo", "growth_gap_hi",
              "sd_gap", "sd_gap_lo", "sd_gap_hi")], row.names = FALSE, digits = 3)
write.csv(ERA, file.path(OUT, "historical_by_era.csv"), row.names = FALSE)
msg("\ngrowth_gap > 0 means high-compliance country-years were followed by faster growth.")
msg("sd_gap < 0 means they were followed by steadier growth.")

# take-off rule, every country with data: mean growth t..t+7 >= 3.5%, exceeding t-7..t-1
# by >= 2pp, and gdp/capita at t+7 above its prior peak (a real take-off, not a rebound)
hdr("BLOCK 2 -- growth take-offs, identified by rule on the whole panel")

FWD8 <- sapply(0:7, function(h) leadv(d$gdp_growth, h))
BWD7 <- sapply(1:7, function(h) lagv(d$gdp_growth, h))
g_post <- ifelse(rowSums(is.na(FWD8)) == 0, rowMeans(FWD8), NA)
g_pre  <- ifelse(rowSums(is.na(BWD7)) == 0, rowMeans(BWD7), NA)
pc_t7  <- leadv(d$gdp_pc, 7)
PCB    <- cbind(d$gdp_pc, sapply(1:7, function(h) lagv(d$gdp_pc, h)))
pc_max_pre <- ifelse(rowSums(is.na(PCB)) == 0, apply(PCB, 1, max), NA)

d$accel <- as.numeric(!is.na(g_post) & !is.na(g_pre) & !is.na(pc_t7) & !is.na(pc_max_pre) &
                      g_post >= 3.5 & (g_post - g_pre) >= 2 & pc_t7 > pc_max_pre)
d$g_post <- g_post; d$g_pre <- g_pre

cand <- !is.na(g_post) & !is.na(g_pre) & !is.na(pc_t7) & !is.na(pc_max_pre)
msg("country-years eligible for the rule: ", sum(cand),
    "   take-off starts identified: ", sum(d$accel[cand] == 1),
    "   countries with at least one: ", length(unique(d$country[d$accel == 1 & cand])))

AC <- d[which(d$accel == 1), c("country", "year", "status", "g_pre", "g_post", "K3")]
AC <- AC[order(AC$country, AC$year), ]
write.csv(AC, file.path(OUT, "accelerations.csv"), row.names = FALSE)
msg("\nfirst take-off year per country (rule applied to all countries):")
first <- AC[!duplicated(AC$country), ]
print(head(first[order(first$year), ], 40), row.names = FALSE, digits = 3)

hdr("BLOCK 3 -- compliance around take-offs against everything else")

# first take-off starts only, so a long acceleration is not counted many times
d$evt <- NA_integer_
fs <- first[, c("country", "year")]
for (i in seq_len(nrow(fs))) {
  s <- d$country == fs$country[i] & d$year >= fs$year[i] - 5 & d$year <= fs$year[i] + 7
  d$evt[s] <- d$year[s] - fs$year[i]
}

boot_mean <- function(x, cid) {
  ok <- !is.na(x); x <- x[ok]; cid <- cid[ok]
  if (!length(x)) return(c(NA, NA, NA, 0))
  cf <- factor(cid); i <- as.integer(cf); C <- nlevels(cf)
  Si <- as.vector(rowsum(x, i, reorder = TRUE)); Ni <- tabulate(i, nbins = C)
  M <- stats::rmultinom(B_BOOT, size = C, prob = rep(1 / C, C))
  bt <- as.vector(crossprod(M, Si)) / as.vector(crossprod(M, Ni))
  q <- quantile(bt, c(.025, .975), na.rm = TRUE)
  c(sum(Si) / sum(Ni), unname(q[1]), unname(q[2]), length(x))
}

rows3 <- list()
for (w in list(c(-5, -1), c(0, 3), c(4, 7))) {
  s <- which(!is.na(d$evt) & d$evt >= w[1] & d$evt <= w[2] & !is.na(d$K3))
  a <- boot_mean(d$K3[s], d$country[s])
  rows3[[length(rows3) + 1]] <- data.frame(
    window = sprintf("t%+d..t%+d", w[1], w[2]), mean_K3 = a[1], ci_lo = a[2], ci_hi = a[3],
    n = a[4], n_countries = length(unique(d$country[s])))
}
s <- which(is.na(d$evt) & !is.na(d$K3))
a <- boot_mean(d$K3[s], d$country[s])
rows3[[length(rows3) + 1]] <- data.frame(window = "all other country-years",
  mean_K3 = a[1], ci_lo = a[2], ci_hi = a[3], n = a[4],
  n_countries = length(unique(d$country[s])))
AC2 <- do.call(rbind, rows3)
print(AC2, row.names = FALSE, digits = 3)
write.csv(AC2, file.path(OUT, "acceleration_compliance.csv"), row.names = FALSE)

msg("\nThe hypothesis predicts mean_K3 in the run-up (t-5..t-1) and during the")
msg("take-off to exceed the panel-wide level, with intervals that do not overlap it.")

msg("\nby indicator, share inside the corridor:")
for (ind in SET3) {
  v <- d[[paste0("in_", ind)]]
  r <- sapply(list(pre = which(!is.na(d$evt) & d$evt < 0), dur = which(!is.na(d$evt) & d$evt >= 0),
                   oth = which(is.na(d$evt))),
              function(s) round(100 * mean(v[s], na.rm = TRUE), 1))
  msg(sprintf("  %-16s run-up %5.1f%%   during %5.1f%%   elsewhere %5.1f%%", ind, r[1], r[2], r[3]))
}

hdr("BLOCK 4 -- the cases in the hypothesis, shown but not used as evidence")

CASES <- c("Korea, Rep.", "China", "Singapore", "Japan", "Germany", "Malaysia",
           "Thailand", "Indonesia", "Chile", "Poland", "Ireland", "Botswana")
tr <- d[d$country %in% CASES & !is.na(d$K3),
        c("country", "year", "status", "gdp_pc", "gdp_growth", "inflation",
          "broad_money_gdp", "trade_balance", "K3")]
write.csv(tr, file.path(OUT, "case_trajectories.csv"), row.names = FALSE)

for (cc in CASES) {
  x <- d[d$country == cc, ]
  if (!nrow(x) || all(is.na(x$K3))) { msg("  ", cc, ": no scored years"); next }
  dec <- cut(x$year, breaks = c(1959, 1969, 1979, 1989, 1999, 2009, 2023),
             labels = c("60s", "70s", "80s", "90s", "00s", "10s+"))
  k  <- round(tapply(x$K3, dec, mean, na.rm = TRUE), 2)
  g  <- round(tapply(x$gdp_growth, dec, mean, na.rm = TRUE), 1)
  msg(sprintf("  %-14s  K3 by decade: %s", cc,
              paste(sprintf("%s=%s", names(k), ifelse(is.na(k), " . ", k)), collapse = "  ")))
  msg(sprintf("  %-14s  growth      : %s", "",
              paste(sprintf("%s=%s", names(g), ifelse(is.na(g), " . ", g)), collapse = "  ")))
}

hdr("BLOCK 5 -- series written out for F2 and F4")

# F2: adherence by year, all countries and a balanced set. the scored sample grows 24 
#to 136 countries, so a rising line could be composition. the balanced set fixes them

sc  <- d[!is.na(d$K3), ]
tab <- table(sc$country); bal <- names(tab)[tab >= 50]
yr  <- sort(unique(sc$year)); ch <- as.character(yr)
AY <- data.frame(year = yr,
  n_all           = as.integer(tapply(sc$K3, sc$year, length)[ch]),
  n_countries_all = as.integer(tapply(sc$country, sc$year, function(v) length(unique(v)))[ch]),
  mean_K3_all     = as.numeric(tapply(sc$K3, sc$year, mean)[ch]),
  pct_all3_all    = 100 * as.numeric(tapply(sc$K3 == 3, sc$year, mean)[ch]))
scb <- sc[sc$country %in% bal, ]
AY$n_balanced       <- as.integer(tapply(scb$K3, scb$year, length)[ch])
AY$mean_K3_balanced <- as.numeric(tapply(scb$K3, scb$year, mean)[ch])
write.csv(AY, file.path(OUT, "adherence_by_year.csv"), row.names = FALSE)
msg("balanced set: ", length(bal), " countries scored in at least 50 of ", length(yr), " years")
msg("mean score, all countries:      ", round(AY$mean_K3_all[1], 2), " (", yr[1], ") -> ",
    round(tail(AY$mean_K3_all, 1), 2), " (", tail(yr, 1), ")")
msg("mean score, balanced countries: ", round(AY$mean_K3_balanced[1], 2), " (", yr[1], ") -> ",
    round(tail(AY$mean_K3_balanced, 1), 2), " (", tail(yr, 1), ")")

# shock detection: the score and each coridor year by year, so each fall is attributed
BY <- do.call(rbind, lapply(yr, function(y) {
  z <- sc[sc$year == y, ]
  data.frame(year = y, n = nrow(z),
             score = mean(z$K3),
             inflation = 100 * mean(z$in_inflation),
             broad_money_gdp = 100 * mean(z$in_broad_money_gdp),
             trade_balance = 100 * mean(z$in_trade_balance))
}))
BY$d_score <- c(NA, diff(BY$score))
for (v in c("inflation", "broad_money_gdp", "trade_balance"))
  BY[[paste0("d_", v)]] <- c(NA, diff(BY[[v]]))
write.csv(BY, file.path(OUT, "adherence_by_year_by_indicator.csv"), row.names = FALSE)
msg("\nten largest year-on-year falls in the score:")
print(head(BY[order(BY$d_score), c("year", "score", "d_score",
      "d_inflation", "d_broad_money_gdp", "d_trade_balance", "n")], 10),
      row.names = FALSE, digits = 3)

# F4 left: mean score in event time around the first take-off
ET <- do.call(rbind, lapply(sort(unique(d$evt[!is.na(d$evt)])), function(e) {
  s <- which(d$evt == e & !is.na(d$K3))
  if (length(s) < 20) return(NULL)
  a <- boot_mean(d$K3[s], d$country[s])
  data.frame(evt = e, mean_K3 = a[1], ci_lo = a[2], ci_hi = a[3], n = a[4],
             n_countries = length(unique(d$country[s])))
}))
s <- which(is.na(d$evt) & !is.na(d$K3)); a <- boot_mean(d$K3[s], d$country[s])
ET$panel_ref <- a[1]; ET$panel_ref_lo <- a[2]; ET$panel_ref_hi <- a[3]
write.csv(ET, file.path(OUT, "event_time_series.csv"), row.names = FALSE)
print(ET, row.names = FALSE, digits = 3)

# F4 right: the named cases by decade
cd <- d[d$country %in% CASES & !is.na(d$K3), ]
cd$decade <- 10 * (cd$year %/% 10)
CD <- aggregate(cbind(K3, gdp_growth) ~ country + decade, cd, mean)
n_dec <- aggregate(K3 ~ country + decade, cd, length); names(n_dec)[3] <- "n_years"
CD <- merge(CD, n_dec)
CD <- CD[order(CD$country, CD$decade), ]
write.csv(CD, file.path(OUT, "case_decade.csv"), row.names = FALSE)
msg("\ncase trajectories by decade written for ", length(unique(CD$country)), " countries")

msg("\ndone: ", format(Sys.time()))
