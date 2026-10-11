#!/usr/bin/env Rscript
# wdi-derived only. cbie, crises, wgi not redistributed

BASE <- here::here()
CORR <- list(
  Developed  = list(inflation = c(1, 3), real_interest = c(0, 2), broad_money_gdp = c(70, 120),
                    trade_balance = c(-3, 6), unemployment = c(3, 8)),
  Developing = list(inflation = c(2, 5), real_interest = c(0, 3), broad_money_gdp = c(40, 80),
                    trade_balance = c(-8, 3), unemployment = c(5, 12)))

flag <- function(D, ind) {
  lo <- ifelse(D$status == "Developed", CORR$Developed[[ind]][1], CORR$Developing[[ind]][1])
  hi <- ifelse(D$status == "Developed", CORR$Developed[[ind]][2], CORR$Developing[[ind]][2])
  as.integer(D[[ind]] >= lo & D[[ind]] < hi) }
lab <- function(s) ifelse(s == "Developing", "Transforming", s)

L <- read.csv(file.path(BASE, "data/wdi_long_1960_2023.csv"))
L$status <- ifelse(is.na(L$gdp_pc), NA, ifelse(L$gdp_pc >= 12696, "Developed", "Developing"))
S3 <- c("inflation", "broad_money_gdp", "trade_balance")
for (j in S3) L[[paste0("in_", j)]] <- flag(L, j)
L$K3 <- ifelse(rowSums(is.na(L[, paste0("in_", S3)])) > 0, NA, rowSums(L[, paste0("in_", S3)]))
L$status <- lab(L$status)

outL <- L[, c("iso3c", "country", "year", "status", "gdp_pc", "gdp_growth", S3, paste0("in_", S3), "K3")]
write.csv(outL, file.path(BASE, "data/adherence_panel_long_1960_2023.csv"), row.names = FALSE)

source(file.path(BASE, "R/panel.R")); P <- load_panel(); P$status <- P$development_status
S5 <- c("inflation", "real_interest", "broad_money_gdp", "trade_balance", "unemployment")
for (j in S5) P[[paste0("in_", j)]] <- flag(P, j)
S4 <- setdiff(S5, "real_interest")

P$K4 <- ifelse(rowSums(is.na(P[, paste0("in_", S4)])) > 0, NA, rowSums(P[, paste0("in_", S4)]))
P$K5 <- ifelse(rowSums(is.na(P[, paste0("in_", S5)])) > 0, NA, rowSums(P[, paste0("in_", S5)]))
P$status <- lab(P$status)
outP <- P[, c("iso2c", "country", "year", "status", "gdp_growth", S5, paste0("in_", S5), "K4", "K5")]
write.csv(outP, file.path(BASE, "data/adherence_panel_detailed_2008_2023.csv"), row.names = FALSE)
cat("long:", nrow(outL), "rows,", sum(!is.na(outL$K3)), "scored |",
    "detailed:", nrow(outP), "rows,", sum(!is.na(outP$K4)), "K4 scored,", sum(!is.na(outP$K5)), "K5 scored\n")

av <- sapply(c("Developed", "Transforming"), function(g) mean(sapply(S5, function(j) {
  m <- P[[paste0("in_", j)]]; s <- P$status == g & !is.na(P$gdp_growth) & !is.na(m); mean(m[s]) })))
cat("check vs paper (41.2 / 30.2):", round(100 * av, 1), "\n")
