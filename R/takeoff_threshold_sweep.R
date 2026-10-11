#!/usr/bin/env Rscript
# grid fixed before the run. paper cell 3.5 / 2 / pc

options(width = 160); set.seed(42); B <- 2000L
BASE <- here::here()
OUT  <- file.path(BASE, "data/foundations"); dir.create(OUT, showWarnings = FALSE)
L <- read.csv(file.path(BASE, "data/wdi_long_1960_2023.csv")); L <- L[order(L$country, L$year), ]
L$status <- ifelse(is.na(L$gdp_pc), NA, ifelse(L$gdp_pc >= 12696, "Developed", "Developing"))
CORR <- list(Developed  = list(inflation = c(1, 3), broad_money_gdp = c(70, 120), trade_balance = c(-3, 6)),
             Developing = list(inflation = c(2, 5), broad_money_gdp = c(40, 80),  trade_balance = c(-8, 3)))

M <- sapply(names(CORR$Developed), function(j) {
  lo <- ifelse(L$status == "Developed", CORR$Developed[[j]][1], CORR$Developing[[j]][1])
  hi <- ifelse(L$status == "Developed", CORR$Developed[[j]][2], CORR$Developing[[j]][2])
  as.numeric(L[[j]] >= lo & L[[j]] < hi) })
L$K3 <- ifelse(rowSums(is.na(M)) > 0, NA, rowSums(M))

key <- paste(L$country, L$year)
ld <- function(v, h) v[match(paste(L$country, L$year + h), key)]
lg <- function(v, h) v[match(paste(L$country, L$year - h), key)]
F8 <- sapply(0:7, function(h) ld(L$gdp_growth, h)); B7 <- sapply(1:7, function(h) lg(L$gdp_growth, h))
gpost <- ifelse(rowSums(is.na(F8)) == 0, rowMeans(F8), NA)
gpre  <- ifelse(rowSums(is.na(B7)) == 0, rowMeans(B7), NA)

pc7 <- ld(L$gdp_pc, 7); PCB <- cbind(L$gdp_pc, sapply(1:7, function(h) lg(L$gdp_pc, h)))
pcmax <- ifelse(rowSums(is.na(PCB)) == 0, apply(PCB, 1, max), NA)

gap_for <- function(G, D, pc) {
  acc <- !is.na(gpost) & !is.na(gpre) & gpost >= G & (gpost - gpre) >= D
  if (pc) acc <- acc & !is.na(pc7) & !is.na(pcmax) & pc7 > pcmax
  A <- L[which(acc), c("country", "year")]; A <- A[!duplicated(A$country), ]
  evt <- rep(NA_integer_, nrow(L))
  for (i in seq_len(nrow(A))) {
    s <- L$country == A$country[i] & L$year >= A$year[i] - 5 & L$year <= A$year[i] + 7
    evt[s] <- L$year[s] - A$year[i] }
  pre <- which(!is.na(evt) & evt < 0 & !is.na(L$K3)); oth <- which(is.na(evt) & !is.na(L$K3))
  est <- mean(L$K3[pre]) - mean(L$K3[oth])
  u <- unique(L$country)
  bt <- replicate(B, { s <- sample(u, replace = TRUE); w <- tabulate(match(L$country, s), length(u))
    w <- as.numeric(table(factor(s, levels = u)))[match(L$country, u)]
    sum(L$K3[pre] * w[pre]) / sum(w[pre]) - sum(L$K3[oth] * w[oth]) / sum(w[oth]) })
  data.frame(G = G, D = D, pc_condition = pc, n_takeoffs = nrow(A), runup_mean = mean(L$K3[pre]),
             other_mean = mean(L$K3[oth]), gap = est, lo = quantile(bt, .025, na.rm = TRUE),
             hi = quantile(bt, .975, na.rm = TRUE), row.names = NULL) }

R <- do.call(rbind, lapply(c(TRUE, FALSE), function(pc) do.call(rbind, lapply(c(2.5, 3, 3.5, 4, 4.5, 5), function(G)
  do.call(rbind, lapply(c(1, 1.5, 2, 2.5, 3), function(D) gap_for(G, D, pc)))))))
write.csv(R, file.path(OUT, "D_takeoff_threshold_sweep.csv"), row.names = FALSE)
print(R, digits = 3, row.names = FALSE)
cat("\ncells:", nrow(R), " gap < 0:", sum(R$gap < 0), " CI below zero:", sum(R$hi < 0),
    " CI above zero:", sum(R$lo > 0), "\n")
cat("paper cell (3.5, 2, pc): "); print(R[R$G == 3.5 & R$D == 2 & R$pc_condition, ], digits = 3, row.names = FALSE)
