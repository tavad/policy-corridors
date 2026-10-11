# sourced by foundations_validation.R, crisis_table.R

CORR <- list(
  Developed  = list(inflation = c(1, 3), real_interest = c(0, 2), broad_money_gdp = c(70, 120),
                    trade_balance = c(-3, 6), unemployment = c(3, 8)),
  Developing = list(inflation = c(2, 5), real_interest = c(0, 3), broad_money_gdp = c(40, 80),
                    trade_balance = c(-8, 3), unemployment = c(5, 12)))

UNIT <- c(inflation = 1, real_interest = 1, broad_money_gdp = 10, trade_balance = 1.5,
          unemployment = 1)
SET3 <- c("inflation", "broad_money_gdp", "trade_balance")
SET4 <- c("inflation", "broad_money_gdp", "trade_balance", "unemployment")
SET5 <- c(SET4, "real_interest")

L <- read.csv(file.path(BASE, "data/wdi_long_1960_2023.csv"))
L <- L[order(L$country, L$year), ]
L$status <- ifelse(is.na(L$gdp_pc), NA, ifelse(L$gdp_pc >= 12696, "Developed", "Developing"))
Lkey <- paste(L$country, L$year)

Llead <- function(v, h) v[match(paste(L$country, L$year + h), Lkey)]
Llag  <- function(v, h) v[match(paste(L$country, L$year - h), Lkey)]
FW5 <- sapply(1:5, function(h) Llead(L$gdp_growth, h))
L$sd_fwd5 <- ifelse(rowSums(is.na(FW5)) == 0, apply(FW5, 1, sd), NA)

FWD8 <- sapply(0:7, function(h) Llead(L$gdp_growth, h))
BWD7 <- sapply(1:7, function(h) Llag(L$gdp_growth, h))
g_post <- ifelse(rowSums(is.na(FWD8)) == 0, rowMeans(FWD8), NA)
g_pre  <- ifelse(rowSums(is.na(BWD7)) == 0, rowMeans(BWD7), NA)

pc7 <- Llead(L$gdp_pc, 7)
PCB <- cbind(L$gdp_pc, sapply(1:7, function(h) Llag(L$gdp_pc, h)))
pcmax <- ifelse(rowSums(is.na(PCB)) == 0, apply(PCB, 1, max), NA)
acc <- !is.na(g_post) & !is.na(g_pre) & !is.na(pc7) & !is.na(pcmax) &
  g_post >= 3.5 & (g_post - g_pre) >= 2 & pc7 > pcmax

AC <- L[which(acc), c("country", "year")]; AC <- AC[!duplicated(AC$country), ]
L$evt <- NA_integer_
for (i in seq_len(nrow(AC))) {
  s <- L$country == AC$country[i] & L$year >= AC$year[i] - 5 & L$year <= AC$year[i] + 7
  L$evt[s] <- L$year[s] - AC$year[i]
}
msg("long panel: ", nrow(L), " rows; take-off first starts: ", nrow(AC))

source(file.path(BASE, "R/panel.R")); P <- load_panel()
P$status <- P$development_status
Pkey <- paste(P$country, P$year)
Plead <- function(v, h) v[match(paste(P$country, P$year + h), Pkey)]

Plag  <- function(v, h) v[match(paste(P$country, P$year - h), Pkey)]
sdw <- function(M) ifelse(rowSums(is.na(M)) == 0, apply(M, 1, sd), NA)
P$sd_f3 <- sdw(sapply(1:3, function(h) Plead(P$gdp_growth, h)))
P$sd_b3 <- sdw(sapply(1:3, function(h) Plag(P$gdp_growth, h)))
iso <- unique(L[, c("iso2c", "iso3c")]); P$iso3c <- iso$iso3c[match(P$iso2c, iso$iso2c)]

inband <- function(x, lo, hi) as.numeric(x >= lo & x < hi)
score <- function(D, C, inds) {
  M <- sapply(inds, function(ind) {
    lo <- ifelse(D$status == "Developed", C$Developed[[ind]][1], C$Developing[[ind]][1])
    hi <- ifelse(D$status == "Developed", C$Developed[[ind]][2], C$Developing[[ind]][2])
    inband(D[[ind]], lo, hi) })
  list(M = M, K = ifelse(rowSums(is.na(M)) > 0, NA_real_, rowSums(M)))
}

twfe_demean <- function(v, a, b, it = 50) {
  for (k in seq_len(it)) { v <- v - ave(v, a); v <- v - ave(v, b) }; v }
TWFE_IT <- 25L
twfe_slope <- function(y, x, cid, yr) {
  ok <- !is.na(y) & !is.na(x); if (sum(ok) < 30) return(NA_real_)
  y <- y[ok]; x <- x[ok]; cid <- cid[ok]; yr <- yr[ok]
  yy <- twfe_demean(y, cid, yr, TWFE_IT); xx <- twfe_demean(x, cid, yr, TWFE_IT)
  sum(xx * yy) / sum(xx * xx) }

auc <- function(score, y) {
  ok <- !is.na(score) & !is.na(y); s <- score[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0); if (n1 < 5 || n0 < 5) return(NA_real_)
  r <- rank(s); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
cboot <- function(fun, cid, B. = B) {
  u <- unique(cid); idx <- split(seq_along(cid), factor(cid, levels = u))
  replicate(B., { s <- sample(length(u), replace = TRUE)
    ii <- unlist(idx[s], use.names = FALSE)
    newid <- rep(seq_along(s), lengths(idx[s]))
    fun(ii, newid) }) }

ci <- function(bt) unname(quantile(bt, c(.025, .975), na.rm = TRUE))
fmt <- function(est, b) sprintf("%.3f [%.3f, %.3f]", est, ci(b)[1], ci(b)[2])

R <- as.data.frame(haven::read_dta(file.path(BASE, "data/ext/cbie_polecon.dta")))
R <- R[order(R$iso_a3, R$year), ]
# s*c2y = crisis in t..t-2: a 1972-73 run predates the file
onset <- function(v, id, yr) { prev <- ave(v, id, FUN = function(z) c(NA, head(z, -1)))
  as.numeric(v == 1 & (is.na(prev) | prev == 0) & yr >= 1974) }

R$on_bank <- onset(R$sbc2y, R$iso_a3, R$year)
R$on_curr <- onset(R$scc2y, R$iso_a3, R$year)
R$on_debt <- onset(R$sdc2y, R$iso_a3, R$year)
# no inflation crises, they are the inflation corridor
R$on_any  <- pmax(R$on_bank, R$on_curr, R$on_debt)
msg("crisis onsets 1974-2017: banking ", sum(R$on_bank), ", currency ", sum(R$on_curr),
    ", debt ", sum(R$on_debt), ", any ", sum(R$on_any),
    " (", length(unique(R$iso_a3)), " economies in file)")

Rk <- paste(R$iso_a3, R$year)
for (v in c("on_any", "on_bank", "on_curr", "on_debt")) {
  L[[v]] <- R[[v]][match(paste(L$iso3c, L$year), Rk)]
}
win <- function(v, hs, lead = TRUE) {
  M <- sapply(hs, function(h) if (lead) Llead(v, h) else Llag(v, h))
  ifelse(rowSums(is.na(M)) > 0, NA, as.numeric(rowSums(M) > 0)) }

L$cr_fwd <- win(L$on_any, 1:3); L$cr_bwd <- win(L$on_any, 1:3, lead = FALSE)
L$cr_fwd_bank <- win(L$on_bank, 1:3)
L$cr_fwd_curr <- win(L$on_curr, 1:3)
L$cr_fwd_debt <- win(L$on_debt, 1:3)
L$inR <- L$iso3c %in% R$iso_a3 & L$year >= 1972 & L$year <= 2017
