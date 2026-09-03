#!/usr/bin/env Rscript

# does joint corridor compliance (corridors met at once) predict forward outcomes?
# K = corridors met at t, complete cases only, so K is not confounded with how many
# series a country reports; gdp_growth is never in K, it is the outcome.
# forward is the test (t -> t+1..t+3); contemporaneous is comparison only.
# SET4 leads; the real lending rate is missing for 44.3% and costs most of the sample.
# every computed test is written out, nothing dropped after its result is seen; BH q-values within each family

options(width = 170, stringsAsFactors = FALSE)
suppressWarnings(RNGversion(getRversion()))

BASE <- here::here()
OUT  <- file.path(BASE, "data")
dir.create(OUT, showWarnings = FALSE)

SEED   <- 42L
B_BOOT <- 20000L
N_PERM <- 2000L
set.seed(SEED)

logcon <- file(file.path(OUT, "joint_run_log.txt"), open = "wt")
sink(logcon, split = TRUE)
on.exit({ sink(); close(logcon) }, add = TRUE)

msg <- function(...) cat(..., "\n", sep = "")
hdr <- function(s) cat("\n", strrep("=", 104), "\n", s, "\n", strrep("=", 104), "\n", sep = "")

hdr("joint_compliance.R -- does simultaneous corridor compliance predict forward outcomes?")
msg("run: ", format(Sys.time()), "   seed=", SEED, "  B_BOOT=", B_BOOT, "  N_PERM=", N_PERM)

hdr("BLOCK 0 -- data, the joint-compliance score, and outcome construction")

source(file.path(BASE, "R/panel.R"))
d <- load_panel()            # 217 economies; the 8 World Bank aggregates dropped
msg("rows=", nrow(d), "  countries=", length(unique(d$country)),
    "  years=", min(d$year), "-", max(d$year))

CORR <- list(
  Developed  = list(inflation = c(1, 3), real_interest = c(0, 2), broad_money_gdp = c(70, 120),
                    trade_balance = c(-3, 6), unemployment = c(3, 8)),
  Developing = list(inflation = c(2, 5), real_interest = c(0, 3), broad_money_gdp = c(40, 80),
                    trade_balance = c(-8, 3), unemployment = c(5, 12)))

SETS <- list(SET4 = c("inflation", "broad_money_gdp", "trade_balance", "unemployment"),
             SET5 = c("inflation", "broad_money_gdp", "trade_balance", "unemployment", "real_interest"))

make_K <- function(z, grp, inds) {
  comp <- sapply(inds, function(ind) {
    b <- CORR[[grp]][[ind]]
    as.numeric(z[[ind]] >= b[1] & z[[ind]] < b[2])
  })
  ifelse(rowSums(is.na(comp)) > 0, NA_real_, rowSums(comp))
}

key  <- paste(d$country, d$year)
LEAD <- sapply(1:5, function(h) d$gdp_growth[match(paste(d$country, d$year + h), key)])
LAGG <- sapply(0:4, function(h) d$gdp_growth[match(paste(d$country, d$year - h), key)])

cmean_g  <- tapply(d$gdp_growth, d$country, mean, na.rm = TRUE)
d$g_cmean <- as.numeric(cmean_g[d$country])

CUT_ABS <- -2
CUT_DEC <- as.numeric(quantile(d$gdp_growth, 0.10, na.rm = TRUE))
msg("large-negative cuts: absolute growth < ", CUT_ABS, " % ; panel 10th pct = ", round(CUT_DEC, 3), " %")

rowStat <- function(M, f) { out <- rep(NA_real_, nrow(M))
  ok <- rowSums(is.na(M)) == 0; if (any(ok)) out[ok] <- apply(M[ok, , drop = FALSE], 1, f); out }

mk_outcomes <- function(M, ref) list(
  meang  = rowStat(M, mean),
  sd     = rowStat(M, sd),
  mad    = { r <- rowMeans(abs(M - ref)); r[rowSums(is.na(M)) > 0] <- NA; r },
  ming   = rowStat(M, min),
  badabs = { r <- as.numeric(rowSums(M < CUT_ABS) > 0); r[rowSums(is.na(M)) > 0] <- NA; r },
  baddec = { r <- as.numeric(rowSums(M < CUT_DEC) > 0); r[rowSums(is.na(M)) > 0] <- NA; r })

FWD3 <- mk_outcomes(LEAD[, 1:3, drop = FALSE], d$g_cmean)
FWD5 <- mk_outcomes(LEAD[, 1:5, drop = FALSE], d$g_cmean)
CON3 <- mk_outcomes(LAGG[, 1:3, drop = FALSE], d$g_cmean)
CON5 <- mk_outcomes(LAGG[, 1:5, drop = FALSE], d$g_cmean)
# pre-trend: K at t against t-3..t-1, strictly before treatment; if K predicts the
# past as well as the future it marks a persistent state and leads nothing
BAK3 <- mk_outcomes(LAGG[, 2:4, drop = FALSE], d$g_cmean)

# theory_sign is the sign one more corridor should have if the framework holds:
# dispersion and tail down, minimum growth up; the mean-growth vs dispersion contrast is the point
OUTCOMES <- data.frame(
  outcome = c("mean_growth_3y", "sd_growth_3y", "sd_growth_5y", "mad_from_country_mean_3y",
              "min_growth_3y", "p_growth_below_-2pct_3y", "p_growth_below_panel_p10_3y"),
  slot    = c("meang", "sd", "sd", "mad", "ming", "badabs", "baddec"),
  window  = c("w3", "w3", "w5", "w3", "w3", "w3", "w3"),
  theory_sign = c(+1, -1, -1, -1, +1, -1, -1))

get_outcome <- function(oi, spec) {
  s <- OUTCOMES$slot[oi]; w <- OUTCOMES$window[oi]
  if (spec == "contemporaneous")  (if (w == "w3") CON3 else CON5)[[s]]
  else if (spec == "backward")    BAK3[[s]]          # 3-year window only
  else                            (if (w == "w3") FWD3 else FWD5)[[s]]
}

pre_covid_ymax <- function(spec, w)
  if (spec == "backward") 2019 else if (w == "w3") 2016 else 2014

# between = pooled OLS slope, within = country fixed-effects slope. the cluster
# bootstrap reduces to per-country multinomial draws, since every term is a sum over
# countries; numerically identical to a resampling loop, fast at 20,000 draws

MULT <- function(C, B) stats::rmultinom(B, size = C, prob = rep(1 / C, C))

boot_p <- function(bt) {
  bt <- bt[!is.na(bt)]; if (!length(bt)) return(NA_real_)
  max(2 * min(mean(bt <= 0), mean(bt >= 0)), 1 / length(bt))
}

fit_between <- function(y, X, cid) {
  cf <- factor(cid); i <- as.integer(cf); C <- nlevels(cf)
  n_i  <- tabulate(i, nbins = C)
  Sx   <- as.vector(rowsum(X,     i, reorder = TRUE))
  Sy   <- as.vector(rowsum(y,     i, reorder = TRUE))
  Sxy  <- as.vector(rowsum(X * y, i, reorder = TRUE))
  Sxx  <- as.vector(rowsum(X * X, i, reorder = TRUE))
  bfun <- function(n, sx, sy, sxy, sxx) {
    den <- sxx - sx^2 / n
    ifelse(n > 1 & den > 0, (sxy - sx * sy / n) / den, NA_real_)
  }
  est <- bfun(sum(n_i), sum(Sx), sum(Sy), sum(Sxy), sum(Sxx))
  M  <- MULT(C, B_BOOT)
  bt <- bfun(as.vector(crossprod(M, n_i)),  as.vector(crossprod(M, Sx)),
             as.vector(crossprod(M, Sy)),   as.vector(crossprod(M, Sxy)),
             as.vector(crossprod(M, Sxx)))
  q <- quantile(bt, c(.025, .975), na.rm = TRUE)
  n <- length(y)
  pn <- replicate(N_PERM, { Xp <- X[sample.int(n)]
    bfun(n, sum(Xp), sum(y), sum(Xp * y), sum(Xp * Xp)) })
  list(est = est, lo = unname(q[1]), hi = unname(q[2]), p = boot_p(bt),
       p_perm = max(mean(abs(pn) >= abs(est), na.rm = TRUE), 1 / N_PERM),
       n_switch = sum(n_i > 0))
}

fit_within <- function(y, X, cid) {
  o <- order(cid); y <- y[o]; X <- X[o]; cid <- cid[o]
  cf <- factor(cid); i <- as.integer(cf); C <- nlevels(cf)
  ni   <- tabulate(i, nbins = C)
  ybar <- as.vector(rowsum(y, i, reorder = TRUE)) / ni
  Xbar <- as.vector(rowsum(X, i, reorder = TRUE)) / ni
  ytil <- y - ybar[i]; Xtil <- X - Xbar[i]
  num_i <- as.vector(rowsum(Xtil * ytil, i, reorder = TRUE))
  den_i <- as.vector(rowsum(Xtil * Xtil, i, reorder = TRUE))
  n_switch <- sum(den_i > 1e-12)              # countries whose own K moves
  if (sum(den_i) <= 0 || n_switch < 5)
    return(list(est = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_,
                p_perm = NA_real_, n_switch = n_switch))
  est <- sum(num_i) / sum(den_i)
  M  <- MULT(C, B_BOOT)
  NU <- as.vector(crossprod(M, num_i)); DE <- as.vector(crossprod(M, den_i))
  bt <- ifelse(DE > 0, NU / DE, NA_real_)
  q  <- quantile(bt, c(.025, .975), na.rm = TRUE)
  n <- length(y)
  pn <- replicate(N_PERM, { Xp <- Xtil[order(i, runif(n))]   # shuffle within country
    sum(Xp * ytil) / sum(Xp * Xp) })
  list(est = est, lo = unname(q[1]), hi = unname(q[2]), p = boot_p(bt),
       p_perm = max(mean(abs(pn) >= abs(est), na.rm = TRUE), 1 / N_PERM), n_switch = n_switch)
}

demean2 <- function(v, ci, yi, tol = 1e-10, maxit = 500) {
  for (it in seq_len(maxit)) {
    v0 <- v
    v <- v - as.vector(rowsum(v, ci, reorder = TRUE) / tabulate(ci))[ci]
    v <- v - as.vector(rowsum(v, yi, reorder = TRUE) / tabulate(yi))[yi]
    if (max(abs(v - v0)) < tol) break
  }
  v
}

MIN_N     <- 60L   # minimum usable observations in a cell
MIN_LEVELS <- 3L   # minimum distinct values of K (dose treatment)

run_cell <- function(y, X, cid, yrv, design, variant) {
  ok <- !is.na(y) & !is.na(X); y <- y[ok]; X <- X[ok]; cid <- cid[ok]; yrv <- yrv[ok]
  base <- data.frame(effect = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p = NA_real_,
                     p_perm = NA_real_, n = length(y), n_treated = sum(X >= 1),
                     mean_K = if (length(X)) mean(X) else NA_real_,
                     n_levels = length(unique(X)),
                     n_countries = length(unique(cid)), n_switchers = NA_integer_, note = "")
  if (base$n < MIN_N || base$n_levels < 2) { base$note <- "insufficient cell"; return(base) }
  if (variant != "yearFE") {
    f <- if (design == "between") fit_between(y, X, cid) else fit_within(y, X, cid)
  } else {
    yi <- as.integer(factor(yrv)); ci <- as.integer(factor(cid))
    if (design == "between") {
      yd <- y - as.vector(rowsum(y, yi, reorder = TRUE) / tabulate(yi))[yi]
      Xd <- X - as.vector(rowsum(X, yi, reorder = TRUE) / tabulate(yi))[yi]
      f <- fit_between(yd, Xd, cid)
    } else {
      f <- fit_within(demean2(y, ci, yi), demean2(X, ci, yi), cid)
    }
  }
  base$effect <- f$est; base$ci_lo <- f$lo; base$ci_hi <- f$hi
  base$p <- f$p; base$p_perm <- f$p_perm; base$n_switchers <- f$n_switch
  if (is.na(f$est)) base$note <- "too few countries with variation in K"
  base
}

hdr("BLOCK 1 -- distribution of the joint-compliance score")

desc <- list()
for (sn in names(SETS)) for (grp in names(CORR)) {
  z <- d[d$development_status == grp, ]
  zi <- which(d$development_status == grp)
  K <- make_K(z, grp, SETS[[sn]])
  for (k in sort(unique(K[!is.na(K)]))) {
    s <- which(K == k)
    desc[[length(desc) + 1]] <- data.frame(
      set = sn, group = grp, K = k, n = length(s),
      n_countries = length(unique(z$country[s])),
      growth_t = mean(z$gdp_growth[s], na.rm = TRUE),
      absdev_t = mean(abs(z$gdp_growth[s] - z$g_cmean[s]), na.rm = TRUE),
      mean_growth_fwd3 = mean(FWD3$meang[zi][s], na.rm = TRUE),
      sd_growth_fwd3   = mean(FWD3$sd[zi][s],    na.rm = TRUE),
      min_growth_fwd3  = mean(FWD3$ming[zi][s],  na.rm = TRUE),
      n_fwd3 = sum(!is.na(FWD3$sd[zi][s])))
  }
  msg("  ", sn, " / ", grp, ": scored country-years = ", sum(!is.na(K)),
      "  countries = ", length(unique(z$country[!is.na(K)])),
      "  countries whose K moves = ",
      sum(tapply(K, z$country, function(v) length(unique(v[!is.na(v)])) > 1), na.rm = TRUE))
}
DESC <- do.call(rbind, desc)
print(DESC, row.names = FALSE, digits = 3)
write.csv(DESC, file.path(OUT, "joint_descriptives.csv"), row.names = FALSE)

hdr("BLOCK 2 -- estimation grid")

TREATS  <- c("dose", "top")        # K linear ; 1{K >= 3}
DESIGNS <- c("between", "within")

variants_for <- function(spec)
  switch(spec,
    forward  = c("none", "yearFE", "preCOVID", "netg", "netg_preCOVID", "sym"),
    backward = c("none", "preCOVID", "sym"),
    "none")

rows <- list(); t0 <- Sys.time()
for (sn in names(SETS)) {
  for (grp in names(CORR)) {
    gsel <- d$development_status == grp
    z <- d[gsel, ]; zi <- which(gsel)
    K <- make_K(z, grp, SETS[[sn]])
    for (spec in c("forward", "contemporaneous", "backward")) {
      for (variant in variants_for(spec)) {
        for (oi in seq_len(nrow(OUTCOMES))) {
          if (spec == "backward" && OUTCOMES$window[oi] == "w5") next   # no 5-year past window
          y <- get_outcome(oi, spec)[zi]
          if (grepl("preCOVID", variant))
            y[z$year > pre_covid_ymax(spec, OUTCOMES$window[oi])] <- NA
          # sym: treatment years whose backward t-3..t-1 and forward t+1..t+3 both lie
          # in 2008-2019, so the pre-trend and the test run on the same years
          if (variant == "sym") y[z$year < 2011 | z$year > 2016] <- NA
          if (grepl("netg", variant)) {
            # netg strips what growth at t already predicts; coefficients are held
            # fixed, which understates uncertainty slightly (flagged in the memo)
            g <- z$gdp_growth; a <- abs(z$gdp_growth - z$g_cmean)
            ok <- !is.na(y) & !is.na(g) & !is.na(a)
            if (sum(ok) > 30) {
              fitn <- lm(y[ok] ~ g[ok] + a[ok])
              yr2 <- rep(NA_real_, length(y)); yr2[ok] <- residuals(fitn); y <- yr2
            } else y[] <- NA
          }
          for (tr in TREATS) {
            X <- if (tr == "dose") K else ifelse(is.na(K), NA, as.numeric(K >= 3))
            for (dg in DESIGNS) {
              r <- run_cell(y, X, z$country, z$year, dg, variant)
              fam <- switch(spec,
                            contemporaneous = "contemporaneous_reference",
                            backward        = paste0("pretrend_", variant),
                            paste0("forward_", variant))
              rows[[length(rows) + 1]] <- cbind(
                data.frame(family = fam, set = sn, group = grp, spec = spec,
                           variant = variant, treatment = tr, design = dg,
                           outcome = OUTCOMES$outcome[oi],
                           theory_sign = OUTCOMES$theory_sign[oi]), r)
            }
          }
        }
      }
    }
    msg("  ... ", sn, " / ", grp, " done (",
        round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s)")
  }
}
res <- do.call(rbind, rows)

hdr("BLOCK 3 -- multiplicity, signs, verdict")

res$supports_theory <- ifelse(is.na(res$effect), NA, sign(res$effect) == res$theory_sign)
res$q_bh <- NA_real_
for (fm in unique(res$family)) {
  s <- res$family == fm & !is.na(res$p)
  res$q_bh[s] <- p.adjust(res$p[s], method = "BH")
  msg("family ", fm, ": ", sum(s), " tests computed, ",
      sum(res$p[s] < .05), " nominal p<.05 (", round(.05 * sum(s), 1), " expected by chance), ",
      sum(res$q_bh[s] < .10, na.rm = TRUE), " with q<.10")
}
res$p_bonf <- ave(res$p, res$family, FUN = function(v) pmin(1, v * sum(!is.na(v))))

msg("\ntotal cells computed: ", nrow(res), "   estimable: ", sum(!is.na(res$effect)))

show <- function(sel, label) {
  s <- res[sel & !is.na(res$effect), ]
  if (!nrow(s)) return(invisible())
  cat("\n---- ", label, " ----\n", sep = "")
  print(s[order(s$group, s$outcome, s$design),
          c("set", "group", "outcome", "treatment", "design", "effect", "ci_lo", "ci_hi",
            "p", "q_bh", "n", "n_countries", "n_switchers", "supports_theory")],
        row.names = FALSE, digits = 3)
}

for (v in c("none", "yearFE", "preCOVID", "netg", "netg_preCOVID", "sym"))
  show(res$family == paste0("forward_", v) & res$set == "SET4" & res$treatment == "dose",
       paste0("SET4, dose (per additional corridor), forward, variant = ", v))
show(res$family == "forward_none" & res$set == "SET4" & res$treatment == "top",
     "SET4, top (K>=3 vs K<=2), forward, no calendar treatment")
show(res$family == "contemporaneous_reference" & res$set == "SET4" & res$treatment == "dose",
     "SET4, dose, CONTEMPORANEOUS -- reference only, not a result (J2)")
for (v in c("none", "preCOVID", "sym"))
  show(res$family == paste0("pretrend_", v) & res$set == "SET4" & res$treatment == "dose",
       paste0("SET4, dose, PRE-TREND (window t-3..t-1), variant = ", v))
show(res$set == "SET5" & res$family == "forward_none" & res$treatment == "dose",
     "SET5 (real rate included), dose, forward, no calendar treatment")

hdr("SIGN TALLY -- forward families only, estimable cells")
for (fm in grep("^forward_", unique(res$family), value = TRUE)) {
  for (dg in DESIGNS) {
    s <- res$family == fm & res$design == dg & !is.na(res$supports_theory)
    msg(sprintf("%-22s %-8s  theory-signed %2d of %2d (%3.0f%%)   nominal p<.05: %2d, of which theory-signed %2d",
                fm, dg, sum(res$supports_theory[s]), sum(s),
                100 * mean(res$supports_theory[s]),
                sum(res$p[s] < .05, na.rm = TRUE),
                sum(res$p[s] < .05 & res$supports_theory[s], na.rm = TRUE)))
  }
}

write.csv(res, file.path(OUT, "joint_results_table.csv"), row.names = FALSE)
msg("\nwritten: ", file.path(OUT, "joint_results_table.csv"))
msg("done: ", format(Sys.time()))
