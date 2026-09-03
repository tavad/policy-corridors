#!/usr/bin/env Rscript
# the corridors as a NECESSARY (not sufficient) condition: does strong performance ever
# occur outside them? necessity predicts an empty cell, and the refutation is a count
# with named counterexamples, not a coefficient.
# success is built from gdp growth only, which is never in the corridor score, so this
# is not the F1 circularity. three nested success definitions, none chosen after the result

options(width = 170, stringsAsFactors = FALSE)
suppressWarnings(RNGversion(getRversion()))

BASE <- here::here()
OUT  <- file.path(BASE, "data/necessity_20260816")
dir.create(OUT, showWarnings = FALSE)

SEED <- 42L; B_BOOT <- 20000L
set.seed(SEED)

logcon <- file(file.path(OUT, "necessity_run_log.txt"), open = "wt")
sink(logcon, split = TRUE); on.exit({ sink(); close(logcon) }, add = TRUE)
msg <- function(...) cat(..., "\n", sep = "")
hdr <- function(s) cat("\n", strrep("=", 104), "\n", s, "\n", strrep("=", 104), "\n", sep = "")

hdr("necessity_test.R -- are the corridors a necessary condition for sustained performance?")
msg("run: ", format(Sys.time()), "   seed=", SEED, "  B_BOOT=", B_BOOT)

source(file.path(BASE, "R/panel.R"))
d <- load_panel()            # 217 economies; the 8 World Bank aggregates dropped

CORR <- list(
  Developed  = list(inflation = c(1, 3), real_interest = c(0, 2), broad_money_gdp = c(70, 120),
                    trade_balance = c(-3, 6), unemployment = c(3, 8)),
  Developing = list(inflation = c(2, 5), real_interest = c(0, 3), broad_money_gdp = c(40, 80),
                    trade_balance = c(-8, 3), unemployment = c(5, 12)))

GROWTH_MIN <- c(Developed = 1, Developing = 5)

SET4 <- c("inflation", "broad_money_gdp", "trade_balance", "unemployment")
SET5 <- c(SET4, "real_interest")

inside <- function(x, b) ifelse(is.na(x), NA, as.numeric(x >= b[1] & x < b[2]))

for (ind in SET5) {
  b <- t(sapply(d$development_status, function(g) CORR[[g]][[ind]]))
  d[[paste0("in_", ind)]] <- ifelse(is.na(d[[ind]]), NA,
                                    as.numeric(d[[ind]] >= b[, 1] & d[[ind]] < b[, 2]))
  d[[paste0("side_", ind)]] <- ifelse(is.na(d[[ind]]), NA,
                                      ifelse(d[[ind]] < b[, 1], "below",
                                             ifelse(d[[ind]] >= b[, 2], "above", "inside")))
}
d$K4 <- ifelse(rowSums(is.na(d[, paste0("in_", SET4)])) > 0, NA,
               rowSums(d[, paste0("in_", SET4)]))

msg("rows=", nrow(d), "  countries=", length(unique(d$country)),
    "  scored (K4) country-years=", sum(!is.na(d$K4)))

hdr("BLOCK 1 -- 5-year spells, three success definitions")

key <- paste(d$country, d$year)
W <- sapply(0:4, function(h) d$gdp_growth[match(paste(d$country, d$year + h), key)])
KW <- sapply(0:4, function(h) d$K4[match(paste(d$country, d$year + h), key)])

ok_w <- rowSums(is.na(W)) == 0
g_mean <- ifelse(ok_w, rowMeans(W), NA)
g_min  <- ifelse(ok_w, apply(W, 1, min), NA)
g_sd   <- ifelse(ok_w, apply(W, 1, sd),  NA)
thr    <- GROWTH_MIN[d$development_status]

sd_med <- tapply(g_sd[ok_w], d$development_status[ok_w], median, na.rm = TRUE)

S1 <- as.numeric(g_mean >= thr)                                        # sustained
S2 <- as.numeric(g_mean >= thr & g_min > 0)                            # + no contraction
S3 <- as.numeric(g_mean >= thr & g_min > 0 &
                 g_sd < sd_med[d$development_status])                  # + steady

# S4 is the framework own success: growth plus a rising export share and rising
# employment; growth alone is what it calls an incomplete indicator
xp0 <- d$exports_goods_services
xp4 <- xp0[match(paste(d$country, d$year + 4), key)]
em0 <- d$employment
em4 <- em0[match(paste(d$country, d$year + 4), key)]
d_exp <- xp4 - xp0
d_emp <- em4 - em0
S4 <- as.numeric(g_mean >= thr & d_exp > 0 & d_emp > 0)                # qualitative
S5 <- as.numeric(g_mean >= thr & d_exp > 0 & d_emp > 0 & g_min > 0)    # + no contraction

SUCCESS <- list(S1_sustained = S1, S2_no_contraction = S2, S3_steady = S3,
                S4_qualitative = S4, S5_qualitative_stable = S5)

# "compliant" = mean score >= 2.5, at least half the corridors held on average
K_mean <- ifelse(rowSums(is.na(KW)) > 0, NA, rowMeans(KW))
COMPLY <- list(K_ge_2.5 = as.numeric(K_mean >= 2.5),
               K_ge_3   = as.numeric(K_mean >= 3),
               K_eq_4   = as.numeric(K_mean >= 3.8))

use <- ok_w & !is.na(K_mean)
msg("complete 5-year spells with a score: ", sum(use),
    "   countries: ", length(unique(d$country[use])))
for (nm in names(SUCCESS))
  msg(sprintf("  base rate of %-18s : Developed %5.1f%%   Developing %5.1f%%", nm,
              100 * mean(SUCCESS[[nm]][use & d$development_status == "Developed"], na.rm = TRUE),
              100 * mean(SUCCESS[[nm]][use & d$development_status == "Developing"], na.rm = TRUE)))
for (nm in names(COMPLY))
  msg(sprintf("  base rate of compliance %-10s: Developed %5.1f%%   Developing %5.1f%%", nm,
              100 * mean(COMPLY[[nm]][use & d$development_status == "Developed"], na.rm = TRUE),
              100 * mean(COMPLY[[nm]][use & d$development_status == "Developing"], na.rm = TRUE)))

# country-cluster bootstrap on a share
boot_share <- function(num, den, cid) {
  ok <- !is.na(num) & !is.na(den); num <- num[ok]; den <- den[ok]; cid <- cid[ok]
  if (!length(num) || sum(den) == 0) return(c(NA, NA, NA, 0, 0))
  cf <- factor(cid); i <- as.integer(cf); C <- nlevels(cf)
  Ni <- as.vector(rowsum(num, i, reorder = TRUE))
  Di <- as.vector(rowsum(den, i, reorder = TRUE))
  est <- sum(Ni) / sum(Di)
  M  <- stats::rmultinom(B_BOOT, size = C, prob = rep(1 / C, C))
  bt <- as.vector(crossprod(M, Ni)) / as.vector(crossprod(M, Di))
  q  <- quantile(bt, c(.025, .975), na.rm = TRUE)
  c(est, unname(q[1]), unname(q[2]), sum(den), sum(num))
}

hdr("BLOCK 2 -- N1: is compliance near-universal among successful spells?")

rows <- list()
for (grp in c("Developed", "Developing")) for (sn in names(SUCCESS)) for (cn in names(COMPLY)) {
  s <- use & d$development_status == grp & !is.na(SUCCESS[[sn]]) & !is.na(COMPLY[[cn]])
  S <- SUCCESS[[sn]][s]; C <- COMPLY[[cn]][s]; cid <- d$country[s]
  a <- boot_share(C * S, S, cid)           # P(compliant | success)
  b <- boot_share(C, rep(1, length(C)), cid)  # P(compliant), base rate
  e <- boot_share(S * (1 - C), rep(1, length(C)), cid)  # share of ALL spells that are successful AND non-compliant
  rows[[length(rows) + 1]] <- data.frame(
    group = grp, success = sn, compliance = cn,
    p_comply_given_success = a[1], ci_lo = a[2], ci_hi = a[3],
    n_success = a[4], n_success_and_compliant = a[5],
    p_comply_base = b[1],
    lift = a[1] - b[1],
    n_counterexamples = e[5],
    pct_counterexamples = 100 * e[1])
}
N1 <- do.call(rbind, rows)
print(N1, row.names = FALSE, digits = 3)
write.csv(N1, file.path(OUT, "necessity_results.csv"), row.names = FALSE)

msg("\nReading: necessity requires p_comply_given_success near 1.  A value close to")
msg("p_comply_base means success carries no information about compliance at all.")

hdr("BLOCK 3 -- N2: per-indicator necessity, country-years")

rows2 <- list()
for (grp in c("Developed", "Developing")) for (ind in SET5) {
  iv <- d[[paste0("in_", ind)]]
  for (sn in names(SUCCESS)) {
    s <- use & d$development_status == grp & !is.na(iv) & !is.na(SUCCESS[[sn]])
    if (sum(s) < 30) next
    S <- SUCCESS[[sn]][s]; I <- iv[s]; cid <- d$country[s]
    a <- boot_share(I * S, S, cid)
    b <- boot_share(I, rep(1, length(I)), cid)
    rows2[[length(rows2) + 1]] <- data.frame(
      group = grp, indicator = ind,
      corridor = sprintf("[%g, %g)", CORR[[grp]][[ind]][1], CORR[[grp]][[ind]][2]),
      success = sn,
      p_inside_given_success = a[1], ci_lo = a[2], ci_hi = a[3],
      p_inside_base = b[1], lift = a[1] - b[1],
      n_success = a[4], n_success_inside = a[5])
  }
}
N2 <- do.call(rbind, rows2)
print(N2[N2$success == "S2_no_contraction",
         c("group", "indicator", "corridor", "p_inside_given_success", "ci_lo", "ci_hi",
           "p_inside_base", "lift", "n_success")], row.names = FALSE, digits = 3)
write.csv(N2, file.path(OUT, "necessity_by_indicator.csv"), row.names = FALSE)

hdr("BLOCK 4 -- N3: is a crossing of the floor the same as a crossing of the ceiling?")

rows3 <- list()
for (grp in c("Developed", "Developing")) for (ind in SET5) {
  sd_v <- d[[paste0("side_", ind)]]
  s0 <- use & d$development_status == grp & !is.na(sd_v)
  for (sd_lab in c("below", "inside", "above")) {
    s <- s0 & sd_v == sd_lab
    if (sum(s) < 20) next
    rows3[[length(rows3) + 1]] <- data.frame(
      group = grp, indicator = ind, side = sd_lab, n = sum(s),
      n_countries = length(unique(d$country[s])),
      mean_growth_5y = mean(g_mean[s], na.rm = TRUE),
      sd_growth_5y   = mean(g_sd[s],   na.rm = TRUE),
      min_growth_5y  = mean(g_min[s],  na.rm = TRUE),
      pct_S2         = 100 * mean(S2[s], na.rm = TRUE))
  }
}
N3 <- do.call(rbind, rows3)
print(N3, row.names = FALSE, digits = 3)
write.csv(N3, file.path(OUT, "asymmetry_by_side.csv"), row.names = FALSE)

hdr("BLOCK 5 -- successful spells that were NOT compliant, named")

s <- use & !is.na(S2) & S2 == 1 & !is.na(K_mean)
CE <- data.frame(country = d$country[s], group = d$development_status[s],
                 spell = paste0(d$year[s], "-", d$year[s] + 4),
                 mean_growth = round(g_mean[s], 2), min_growth = round(g_min[s], 2),
                 sd_growth = round(g_sd[s], 2), mean_K = round(K_mean[s], 2))
CE <- CE[order(CE$mean_K, -CE$mean_growth), ]
write.csv(CE, file.path(OUT, "necessity_counterexamples.csv"), row.names = FALSE)

msg("successful spells (S2) with a score: ", nrow(CE))
for (grp in c("Developed", "Developing")) {
  x <- CE[CE$group == grp, ]
  if (!nrow(x)) next
  msg("\n", grp, ": ", nrow(x), " successful spells, of which ",
      sum(x$mean_K < 2), " ran at mean K < 2 and ", sum(x$mean_K < 1), " at mean K < 1.")
  msg("  least compliant successful spells:")
  print(head(x, 15), row.names = FALSE)
}

msg("\ndone: ", format(Sys.time()))
