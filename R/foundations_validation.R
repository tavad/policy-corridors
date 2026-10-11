#!/usr/bin/env Rscript
options(width = 170, stringsAsFactors = FALSE, warn = 1)
suppressWarnings(RNGversion(getRversion()))
BASE <- here::here()
OUT  <- file.path(BASE, "data/foundations"); dir.create(OUT, showWarnings = FALSE)
SEED <- 42L; B <- as.integer(Sys.getenv("B_BOOT", "2000")); NDRAW <- as.integer(Sys.getenv("NDRAW", "1000"))
set.seed(SEED)

logcon <- file(file.path(OUT, "run_log.txt"), open = "wt")
sink(logcon, split = TRUE)
msg <- function(...) cat(..., "\n", sep = "")
hdr <- function(s) cat("\n", strrep("=", 100), "\n", s, "\n", strrep("=", 100), "\n", sep = "")
hdr("foundations_validation.R"); msg("run ", format(Sys.time()), "  seed ", SEED, "  B ", B, "  draws ", NDRAW)

source(file.path(BASE, "R/external.R"))

findings <- function(C, quick = FALSE) {
  out <- list()
  sc <- score(P, C, SET5)
  okg <- !is.na(P$gdp_growth)
  av <- sapply(c("Developed", "Developing"), function(g) mean(sapply(SET5, function(ind) {
    m <- sc$M[, ind]; s <- P$status == g & okg & !is.na(m); mean(m[s]) })))
  out$adh_dev <- 100 * av[1]; out$adh_tr <- 100 * av[2]
  out$F1 <- av[1] > av[2] && all(av < .5)
  k3 <- score(L, C, SET3)$K
  dec <- cut(L$year, c(1959, 1969, 1979, 1989, 1999, 2009, 2023))
  dm <- tapply(k3, dec, mean, na.rm = TRUE)
  ym <- tapply(k3, L$year, mean, na.rm = TRUE); yy <- as.numeric(names(ym))
  sl80 <- coef(lm(ym[yy >= 1980] ~ yy[yy >= 1980]))[2] * 10
  out$trough <- which.min(dm); out$slope80 <- sl80
  out$F2 <- out$trough %in% c(2, 3) && dm[1] > dm[3] && sl80 > 0
  pre <- mean(k3[!is.na(L$evt) & L$evt < 0], na.rm = TRUE); oth <- mean(k3[is.na(L$evt)], na.rm = TRUE)
  out$takeoff_gap <- pre - oth; out$F3 <- pre < oth
  gap <- sapply(c("Developed", "Developing"), function(g) {
    s <- L$status == g & !is.na(k3) & !is.na(L$sd_fwd5)
    mean(L$sd_fwd5[s & k3 >= 2]) - mean(L$sd_fwd5[s & k3 <= 1]) })
  out$stab_gap <- mean(gap); out$F5 <- out$stab_gap < 0
  s <- L$inR & L$year <= 2014
  out$auc <- auc(-k3[s], L$cr_fwd[s])
  if (!quick) {
    # 2011-16: both windows inside 2008-19
    k4 <- score(P, C, SET4)$K
    s <- P$status == "Developed" & P$year >= 2011 & P$year <= 2016
    out$b_back <- twfe_slope(P$sd_b3[s], k4[s], P$country[s], P$year[s])
    out$b_fwd  <- twfe_slope(P$sd_f3[s], k4[s], P$country[s], P$year[s])
    out$F4 <- isTRUE(out$b_back < out$b_fwd)
  }
  out
}

hdr("BASELINE -- the framework's own corridors")
base <- findings(CORR)
str(base)

hdr("B1 -- bound sensitivity: 1,000 random perturbations + named variants")
perturb <- function() {
  C <- CORR
  for (g in names(C)) for (ind in names(C[[g]])) repeat {
    b <- C[[g]][[ind]] + runif(2, -1, 1) * UNIT[ind]
    if (b[2] - b[1] >= .5 * UNIT[ind]) { C[[g]][[ind]] <- b; break } }
  C }

t0 <- Sys.time()
PR <- do.call(rbind, lapply(seq_len(NDRAW), function(r) { f <- findings(perturb())
  data.frame(draw = r, adh_dev = f$adh_dev, adh_tr = f$adh_tr, F1 = f$F1, trough = f$trough,
             slope80 = f$slope80, F2 = f$F2, takeoff_gap = f$takeoff_gap, F3 = f$F3,
             b_back = f$b_back, b_fwd = f$b_fwd, F4 = f$F4, stab_gap = f$stab_gap, F5 = f$F5,
             auc = f$auc) }))
msg("perturbation run took ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

write.csv(PR, file.path(OUT, "B1_perturbations.csv"), row.names = FALSE)
msg("share of 1,000 perturbed systems in which each finding holds:")
print(round(colMeans(PR[, c("F1", "F2", "F3", "F4", "F5")]), 3))
msg("ranges (5th-95th pct): adh_dev ", paste(round(quantile(PR$adh_dev, c(.05, .95)), 1), collapse = "-"),
    "  adh_tr ", paste(round(quantile(PR$adh_tr, c(.05, .95)), 1), collapse = "-"),
    "  auc ", paste(round(quantile(PR$auc, c(.05, .95), na.rm = TRUE), 3), collapse = "-"))

V <- list()
V$baseline <- CORR
V$stated_dev_inflation_0_2 <- within(CORR, Developed$inflation <- c(0, 2))
V$all_wider_1unit <- CORR; V$all_narrower_halfunit <- CORR

for (g in names(CORR)) for (ind in names(CORR[[g]])) {
  V$all_wider_1unit[[g]][[ind]] <- CORR[[g]][[ind]] + c(-1, 1) * UNIT[ind]
  V$all_narrower_halfunit[[g]][[ind]] <- CORR[[g]][[ind]] + c(.5, -.5) * UNIT[ind] }
# MIP is current account, trade balance as proxy
V$benchmark_anchored <- CORR
V$benchmark_anchored$Developed$inflation  <- c(1, 3)
V$benchmark_anchored$Developing$inflation <- c(1.5, 4.5)
V$benchmark_anchored$Developed$trade_balance <- c(-4, 6)
V$benchmark_anchored$Developed$unemployment  <- c(3, 10)

VV <- do.call(rbind, lapply(names(V), function(n) { f <- findings(V[[n]])
  data.frame(variant = n, adh_dev = f$adh_dev, adh_tr = f$adh_tr, F1 = f$F1, trough = f$trough,
             slope80 = f$slope80, F2 = f$F2, takeoff_gap = f$takeoff_gap, F3 = f$F3,
             b_back = f$b_back, b_fwd = f$b_fwd, F4 = f$F4, stab_gap = f$stab_gap, F5 = f$F5,
             auc = f$auc) }))

print(VV, row.names = FALSE, digits = 3)
write.csv(VV, file.path(OUT, "B1_named_variants.csv"), row.names = FALSE)

hdr("B2 -- falsification test: the framework's bounds against random corridors of equal width")
# same widths, centres in p10-p90
Q <- lapply(c(Developed = "Developed", Developing = "Developing"), function(g)
  sapply(SET5, function(ind) { x <- c(L[[ind]][L$status %in% g], P[[ind]][P$status %in% g])
    quantile(x, c(.1, .9), na.rm = TRUE) }))

randsys <- function() { C <- CORR
  for (g in names(C)) for (ind in names(C[[g]])) { w <- diff(CORR[[g]][[ind]])
    c0 <- runif(1, Q[[g]][1, ind], Q[[g]][2, ind]); C[[g]][[ind]] <- c0 + c(-w, w) / 2 }
  C }
RS <- do.call(rbind, lapply(seq_len(NDRAW), function(r) { f <- findings(randsys(), quick = TRUE)
  data.frame(draw = r, stab_gap = f$stab_gap, auc = f$auc, adh_dev = f$adh_dev, adh_tr = f$adh_tr,
             takeoff_gap = f$takeoff_gap) }))

write.csv(RS, file.path(OUT, "B2_random_systems.csv"), row.names = FALSE)
msg(sprintf("stability gap: framework %.3f; random median %.3f; share of random systems with a MORE negative gap: %.3f",
            base$stab_gap, median(RS$stab_gap), mean(RS$stab_gap < base$stab_gap)))
msg(sprintf("crisis AUC:    framework %.3f; random median %.3f; share of random systems with a HIGHER AUC: %.3f",
            base$auc, median(RS$auc, na.rm = TRUE), mean(RS$auc > base$auc, na.rm = TRUE)))
msg(sprintf("adherence dev/tr: framework %.1f/%.1f; random medians %.1f/%.1f",
            base$adh_dev, base$adh_tr, median(RS$adh_dev), median(RS$adh_tr)))

hdr("B3 -- validation: does the score discriminate crisis onsets?")
L$K3 <- score(L, CORR, SET3)$K
for (j in SET3) L[[paste0("in_", j)]] <- score(L, CORR, j)$M[, 1]
s <- which(L$inR & L$year <= 2014 & !is.na(L$K3) & !is.na(L$cr_fwd))

S <- L[s, ]
msg("sample: ", nrow(S), " country-years, ", length(unique(S$country)), " economies, ",
    sum(S$cr_fwd), " with an onset in t+1..t+3 (base rate ", round(mean(S$cr_fwd), 3), ")")
rateK <- tapply(S$cr_fwd, S$K3, mean); nK <- table(S$K3)
msg("onset rate in t+1..t+3 by K3: ", paste(sprintf("K=%s %.3f (n=%d)", names(rateK), rateK, nK), collapse = "; "))

for (g in c("Developed", "Developing")) { Sg <- S[S$status == g, ]
  r <- tapply(Sg$cr_fwd, Sg$K3, mean); n <- table(Sg$K3)
  msg("  ", g, ": ", paste(sprintf("K=%s %.3f (n=%d)", names(r), r, n), collapse = "; ")) }
rd <- function(ii) { y <- S$cr_fwd[ii]; k <- S$K3[ii]; mean(y[k == 3]) - mean(y[k <= 1]) }
e_rd <- rd(seq_len(nrow(S))); bt_rd <- cboot(function(ii, id) rd(ii), S$country)
msg("onset rate, K3 = 3 minus K3 <= 1: ", fmt(e_rd, bt_rd))

B3 <- list()
add3 <- function(label, est, bt, n) B3[[length(B3) + 1]] <<- data.frame(test = label, est = est,
  lo = ci(bt)[1], hi = ci(bt)[2], p_two = 2 * min(mean(bt <= 0, na.rm = TRUE), mean(bt >= 0, na.rm = TRUE)), n = n)
for (o in c("cr_fwd", "cr_fwd_bank", "cr_fwd_curr", "cr_fwd_debt")) {
  e <- auc(-S$K3, S[[o]]); bt <- cboot(function(ii, id) auc(-S$K3[ii], S[[o]][ii]), S$country)
  add3(paste0("AUC -K3 -> ", o), e, bt - .5 + .5 * 0, nrow(S))
  B3[[length(B3)]]$p_two <- 2 * min(mean(bt <= .5), mean(bt >= .5))
  msg(sprintf("AUC(low score -> %-12s) %s", o, fmt(e, bt)))
}

for (j in SET3) { v <- S[[paste0("in_", j)]]
  e <- auc(-v, S$cr_fwd); bt <- cboot(function(ii, id) auc(-v[ii], S$cr_fwd[ii]), S$country)
  add3(paste0("AUC outside ", j, " -> cr_fwd"), e, bt, nrow(S))
  B3[[length(B3)]]$p_two <- 2 * min(mean(bt <= .5), mean(bt >= .5))
  msg(sprintf("AUC(outside %-15s -> any onset) %s", j, fmt(e, bt))) }

sb <- S[!is.na(S$cr_bwd), ]
for (o in c("cr_fwd", "cr_bwd")) {
  e <- twfe_slope(sb[[o]], sb$K3, sb$country, sb$year)
  bt <- cboot(function(ii, id) twfe_slope(sb[[o]][ii], sb$K3[ii], id, sb$year[ii]), sb$country)
  add3(paste0("TWFE LPM ", o, " ~ K3"), e, bt, nrow(sb))
  msg(sprintf("two-way FE, P(onset %s) per extra corridor: %s  (n=%d)",
              ifelse(o == "cr_fwd", "in t+1..t+3", "in t-3..t-1"), fmt(e, bt), nrow(sb)))
}

cm <- aggregate(cbind(K3, cr_fwd) ~ country, data = S, FUN = mean)
e <- cor(cm$K3, cm$cr_fwd, method = "spearman")
bt <- replicate(B, { i <- sample(nrow(cm), replace = TRUE); cor(cm$K3[i], cm$cr_fwd[i], method = "spearman") })
add3("rate diff K3=3 minus K3<=1", e_rd, bt_rd, nrow(S))
add3("between Spearman meanK3 vs crisis freq", e, bt, nrow(cm))
msg("between-economy Spearman, mean K3 vs onset frequency: ", fmt(e, bt))

for (g in c("Developed", "Developing")) { Sg <- S[S$status == g, ]
  e <- auc(-Sg$K3, Sg$cr_fwd); bt <- cboot(function(ii, id) auc(-Sg$K3[ii], Sg$cr_fwd[ii]), Sg$country)
  add3(paste0("AUC -K3 -> cr_fwd, ", g), e, bt, nrow(Sg))
  B3[[length(B3)]]$p_two <- 2 * min(mean(bt <= .5), mean(bt >= .5))
  msg(sprintf("AUC %-10s %s  (n=%d, onsets %d)", g, fmt(e, bt), nrow(Sg), sum(Sg$cr_fwd))) }

B3 <- do.call(rbind, B3); write.csv(B3, file.path(OUT, "B3_crisis_validation.csv"), row.names = FALSE)

hdr("B4 -- institutions: governance (WGI), central-bank independence (CBIE), democracy")
W <- read.csv(file.path(BASE, "data/ext/wgi_1996_2023.csv"))
W <- W[!is.na(W$iso3c) & nchar(W$iso3c) == 3, ]; W <- W[!duplicated(W[, c("iso3c", "year", "indicator")]), ]
Ww <- reshape(W, idvar = c("iso3c", "year"), timevar = "indicator", direction = "wide")
names(Ww) <- sub("value\\.", "", names(Ww))
# GE RQ RL CC only, §3.1
Ww$gov <- rowMeans(Ww[, c("GE", "RQ", "RL", "CC")])
Wk <- paste(Ww$iso3c, Ww$year)

L$gov_l1 <- Ww$gov[match(paste(L$iso3c, L$year - 1), Wk)]
L$gov    <- Ww$gov[match(paste(L$iso3c, L$year), Wk)]
P$gov_l1 <- Ww$gov[match(paste(P$iso3c, P$year - 1), Wk)]
P$K4 <- score(P, CORR, SET4)$K
L$cbie    <- R$CBIE[match(paste(L$iso3c, L$year), Rk)]
L$cbie_l1 <- R$CBIE[match(paste(L$iso3c, L$year - 1), Rk)]
L$democ   <- R$democracy[match(paste(L$iso3c, L$year), Rk)]

B4 <- list()
add4 <- function(label, est, bt, n, nc) B4[[length(B4) + 1]] <<- data.frame(test = label, est = est,
  lo = ci(bt)[1], hi = ci(bt)[2], p_two = 2 * min(mean(bt <= 0, na.rm = TRUE), mean(bt >= 0, na.rm = TRUE)),
  n = n, n_econ = nc)

betw <- function(D, kv, iv, label, by = NULL) {
  for (g in c("All", "Developed", "Developing")) {
    Dg <- if (g == "All") D else D[D$status %in% g, ]
    Dg <- Dg[!is.na(Dg[[kv]]) & !is.na(Dg[[iv]]), ]
    cm <- aggregate(Dg[, c(kv, iv)], list(c = Dg$country), mean); if (nrow(cm) < 10) next
    e <- cor(cm[[kv]], cm[[iv]], method = "spearman")
    bt <- replicate(B, { i <- sample(nrow(cm), replace = TRUE); cor(cm[[kv]][i], cm[[iv]][i], method = "spearman") })
    add4(paste0("between Spearman ", kv, " vs ", iv, " [", label, ", ", g, "]"), e, bt, nrow(Dg), nrow(cm))
    msg(sprintf("between  %-4s vs %-6s %-24s %-10s rho %s  (%d economies)", kv, iv, label, g, fmt(e, bt), nrow(cm))) } }

betw(L[L$year >= 1996, ], "K3", "gov", "long 1996-2023")
betw(transform(P, gov = gov_l1), "K4", "gov", "detailed 2008-2023")
betw(L[L$year <= 2017, ], "K3", "cbie", "long 1972-2017")
betw(L[L$year <= 2017, ], "K3", "democ", "long 1972-2017")

withn <- function(D, kv, xv, label) for (g in c("All", "Developed", "Developing")) {
  Dg <- if (g == "All") D else D[D$status %in% g, ]
  Dg <- Dg[!is.na(Dg[[kv]]) & !is.na(Dg[[xv]]), ]
  e <- twfe_slope(Dg[[kv]], Dg[[xv]], Dg$country, Dg$year)
  bt <- cboot(function(ii, id) twfe_slope(Dg[[kv]][ii], Dg[[xv]][ii], id, Dg$year[ii]), Dg$country)
  add4(paste0("TWFE ", kv, " ~ ", xv, " [", label, ", ", g, "]"), e, bt, nrow(Dg), length(unique(Dg$country)))
  msg(sprintf("within   %-4s ~ %-7s %-24s %-10s b %s  (n=%d, %d economies)",
              kv, xv, label, g, fmt(e, bt), nrow(Dg), length(unique(Dg$country)))) }

withn(L, "K3", "gov_l1", "long 1997-2023")
withn(P, "K4", "gov_l1", "detailed 2008-2023")
L$K3_l1 <- Llag(L$K3, 1); L$gov_f0 <- L$gov
Z <- L[!is.na(L$K3) & !is.na(L$gov_l1) & !is.na(L$K3_l1) & !is.na(L$gov), ]
Z$zK <- (Z$K3 - mean(Z$K3)) / sd(Z$K3); Z$zKl <- (Z$K3_l1 - mean(Z$K3)) / sd(Z$K3)
Z$zG <- (Z$gov - mean(Z$gov)) / sd(Z$gov); Z$zGl <- (Z$gov_l1 - mean(Z$gov)) / sd(Z$gov)

for (g in c("All", "Developing")) { Zg <- if (g == "All") Z else Z[Z$status %in% g, ]
  for (dir in c("gov(t-1) -> K3(t)", "K3(t-1) -> gov(t)")) {
    yv <- if (dir == "gov(t-1) -> K3(t)") "zK" else "zG"; xv <- if (dir == "gov(t-1) -> K3(t)") "zGl" else "zKl"
    e <- twfe_slope(Zg[[yv]], Zg[[xv]], Zg$country, Zg$year)
    bt <- cboot(function(ii, id) twfe_slope(Zg[[yv]][ii], Zg[[xv]][ii], id, Zg$year[ii]), Zg$country)
    add4(paste0("ordering, standardized TWFE ", dir, " [", g, "]"), e, bt, nrow(Zg), length(unique(Zg$country)))
    msg(sprintf("ordering  %-18s %-10s %s  (n=%d)", dir, g, fmt(e, bt), nrow(Zg))) } }

withn(L, "K3", "cbie_l1", "long 1973-2018")
withn(L, "K3", "democ", "long 1972-2017")

R$dC <- R$CBIE_Variation
ev <- do.call(rbind, lapply(split(R, R$iso_a3), function(z) {
  z <- z[!is.na(z$dC) & z$dC >= .10 & z$year >= 1975 & z$year <= 2013, ]
  if (!nrow(z)) NULL else z[which.max(z$dC), c("iso_a3", "year", "dC")] }))
names(ev) <- c("iso3c", "e", "dC")
anypos <- tapply(seq_len(nrow(R)), R$iso_a3, function(ii) R$year[ii][which(R$dC[ii] > 0)])
pre_on <- sapply(seq_len(nrow(ev)), function(i) {
  z <- R[R$iso_a3 == ev$iso3c[i] & R$year >= ev$e[i] - 3 & R$year <= ev$e[i], ]
  any(z$on_any == 1) })

ev$crisis_before <- pre_on
ev$status <- L$status[match(paste(ev$iso3c, ev$e), paste(L$iso3c, L$year))]
msg("reform events (largest dCBIE >= 0.10 per economy, 1975-2013): ", nrow(ev),
    "; preceded by a crisis onset in [e-3,e]: ", sum(ev$crisis_before),
    "; status at event: ", paste(names(table(ev$status)), table(ev$status), collapse = " / "))

Lc <- L[!is.na(L$K3), c("iso3c", "year", "K3")]
K3m <- tapply(Lc$K3, list(Lc$iso3c, Lc$year), mean)
ctrl_ok <- function(c, e) { y <- anypos[[c]]; is.null(y) || !any(y >= e - 5 & y <= e + 10) }
allc <- rownames(K3m); KK <- -5:10
path <- function(i) { c <- ev$iso3c[i]; e <- ev$e[i]
  if (!c %in% allc) return(rep(NA, length(KK)))
  ctr <- allc[sapply(allc, function(cc) cc != c && (!cc %in% names(anypos) || ctrl_ok(cc, e)))]
  sapply(KK, function(k) { yr <- as.character(e + k); if (!yr %in% colnames(K3m)) return(NA)
    K3m[c, yr] - mean(K3m[ctr, yr], na.rm = TRUE) }) }

PM <- t(sapply(seq_len(nrow(ev)), path)); colnames(PM) <- KK
norm <- PM - rowMeans(PM[, as.character(-5:-1)], na.rm = TRUE)
keep <- rowSums(!is.na(PM[, as.character(-5:-1)])) >= 3 & rowSums(!is.na(PM[, as.character(0:10)])) >= 5
msg("events with >=3 pre and >=5 post scored years: ", sum(keep))
es <- function(rows) { N <- norm[rows, , drop = FALSE]
  est <- colMeans(N, na.rm = TRUE)
  bt <- replicate(B, colMeans(N[sample(nrow(N), replace = TRUE), , drop = FALSE], na.rm = TRUE))
  post <- function(M) mean(M[, as.character(0:10)], na.rm = TRUE)
  postb <- replicate(B, { i <- sample(nrow(N), replace = TRUE); post(N[i, , drop = FALSE]) })
  list(path = data.frame(k = KK, est = est, lo = apply(bt, 1, quantile, .025, na.rm = TRUE),
                         hi = apply(bt, 1, quantile, .975, na.rm = TRUE)),
       post = c(post(N), ci(postb)), n = nrow(N)) }

ESall <- es(which(keep)); ESnc <- es(which(keep & !ev$crisis_before))
ESdev <- es(which(keep & ev$status %in% "Developing"))
for (nm in c("ESall", "ESnc", "ESdev")) { x <- get(nm)
  msg(sprintf("%-6s events %3d  mean shift k=0..10 vs k=-5..-1: %.3f [%.3f, %.3f]", nm, x$n, x$post[1], x$post[2], x$post[3]))
  add4(paste0("CBI reform event study, post-pre mean K3 shift [", nm, "]"), x$post[1], c(x$post[2], x$post[3]), x$n, x$n) }
msg("event-time path, all events:"); print(ESall$path, row.names = FALSE, digits = 3)

write.csv(rbind(cbind(set = "all", ESall$path), cbind(set = "no_prior_crisis", ESnc$path),
                cbind(set = "transforming_at_event", ESdev$path)),
          file.path(OUT, "B4_cbi_event_study.csv"), row.names = FALSE)
write.csv(ev, file.path(OUT, "B4_cbi_events.csv"), row.names = FALSE)
B4 <- do.call(rbind, B4); B4$p_two[grepl("event study", B4$test)] <- NA

write.csv(B4, file.path(OUT, "B4_institutions.csv"), row.names = FALSE)
hdr("done"); msg(format(Sys.time()))
sink(); close(logcon)
