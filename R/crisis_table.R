#!/usr/bin/env Rscript
options(width = 170, stringsAsFactors = FALSE, warn = 1)
suppressWarnings(RNGversion(getRversion()))
BASE <- here::here()
OUT  <- file.path(BASE, "data/foundations"); dir.create(OUT, showWarnings = FALSE)
SEED <- 42L; B <- as.integer(Sys.getenv("B_BOOT", "2000"))
set.seed(SEED)

logcon <- file(file.path(OUT, "crisis_table_log.txt"), open = "wt")
sink(logcon, split = TRUE)
msg <- function(...) cat(..., "\n", sep = "")
hdr <- function(s) cat("\n", strrep("=", 100), "\n", s, "\n", strrep("=", 100), "\n", sep = "")
hdr("crisis_table.R"); msg("run ", format(Sys.time()), "  seed ", SEED, "  B ", B)

source(file.path(BASE, "R/external.R"))

L$K3 <- score(L, CORR, SET3)$K
L$K2 <- score(L, CORR, c("broad_money_gdp", "trade_balance"))$K
L$in_inf <- score(L, CORR, "inflation")$M[, 1]
S <- L[which(L$inR & L$year <= 2014 & !is.na(L$K3) & !is.na(L$cr_fwd)), ]

S$yr_rate <- ave(S$cr_fwd, S$year)
cellrate <- function(D, kv, ov = "cr_fwd", adj = FALSE) do.call(rbind, lapply(sort(unique(D[[kv]])), function(k) {
  Dk <- D[D[[kv]] == k & !is.na(D[[ov]]), ]; y <- if (adj) Dk[[ov]] - Dk$yr_rate else Dk[[ov]]
  bt <- cboot(function(ii, id) mean(y[ii]), Dk$country)
  data.frame(k = k, rate = mean(y), lo = ci(bt)[1], hi = ci(bt)[2], n = nrow(Dk),
             economies = length(unique(Dk$country)), onsets = sum(Dk[[ov]])) }))

OUTT <- list()
put <- function(label, x) { x <- cbind(table = label, x); OUTT[[length(OUTT) + 1]] <<- x
  msg("\n", label); print(x[, -1], row.names = FALSE, digits = 3) }
for (g in c("All", "Developed", "Developing")) { D <- if (g == "All") S else S[S$status == g, ]
  put(paste0("T1 K3 -> onset t+1..t+3, ", g), cellrate(D, "K3"))
  put(paste0("T2 K2 (no inflation) -> onset, ", g), cellrate(D, "K2"))
  put(paste0("T2b inflation corridor alone -> onset, ", g), cellrate(D, "in_inf")) }

Dd <- S[S$status == "Developing", ]
for (o in c("cr_fwd_bank", "cr_fwd_curr", "cr_fwd_debt")) put(paste0("T3 K3 -> ", o, ", Developing"), cellrate(Dd, "K3", o))
put("T4 K3 -> onset, Developing, 1974-1989", cellrate(Dd[Dd$year <= 1989, ], "K3"))
put("T4 K3 -> onset, Developing, 1990-2014", cellrate(Dd[Dd$year >= 1990, ], "K3"))
put("T5 K3 -> onset minus same-year rate, Developing", cellrate(Dd, "K3", adj = TRUE))
put("T6 K3 at t <- onset in t-3..t-1 (mirror), Developing", cellrate(Dd, "K3", "cr_bwd"))
put("T6 K3 at t <- onset in t-3..t-1 (mirror), All", cellrate(S, "K3", "cr_bwd"))

k3 <- Dd[Dd$K3 == 3, ]; tb <- sort(table(k3$country), decreasing = TRUE)
msg("\nT7 developing K3=3 country-years: ", nrow(k3), " from ", length(tb), " economies; top 10:")
print(head(tb, 10)); msg("share of K3=3 years from top 5 economies: ", round(sum(head(tb, 5)) / nrow(k3), 3))
msg("onsets among developing K3=3 years: ", sum(k3$cr_fwd), " in ",
    paste(unique(k3$country[k3$cr_fwd == 1]), collapse = ", "))
# K3=3 cell: 4 onsets, 2 economies
loo <- sapply(names(tb), function(c) mean(k3$cr_fwd[k3$country != c]))
msg("leave-one-economy-out range of the developing K3=3 rate: ", paste(round(range(loo), 3), collapse = " - "))

write.csv(do.call(rbind, OUTT), file.path(OUT, "C_crisis_table.csv"), row.names = FALSE)
sink(); close(logcon)
