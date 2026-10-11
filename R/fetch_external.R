#!/usr/bin/env Rscript
suppressMessages(library(jsonlite))
BASE <- here::here()
EXT  <- file.path(BASE, "data/ext")
dir.create(EXT, showWarnings = FALSE, recursive = TRUE)

codes <- c("GE", "RQ", "RL", "CC", "VA", "PV")
out <- list()
for (k in codes) {
  u <- sprintf(paste0("https://api.worldbank.org/v2/country/all/indicator/GOV_WGI_%s.EST",
                      "?source=3&format=json&per_page=20000&date=1996:2023"), k)
  j <- fromJSON(u)
  x <- j[[2]]
  out[[k]] <- data.frame(iso3c = x$countryiso3code, year = as.integer(x$date),
                         indicator = k, value = x$value)
  cat(k, nrow(out[[k]]), sum(!is.na(x$value)), "\n")
}

w <- do.call(rbind, out)
write.csv(w, file.path(EXT, "wgi_1996_2023.csv"), row.names = FALSE)

# format=original or dataverse serves its tab export
download.file("https://dataverse.harvard.edu/api/access/datafile/6062520?format=original",
              file.path(EXT, "cbie_polecon.dta"), mode = "wb", quiet = TRUE)

writeLines(paste("WGI fetched", format(Sys.time()), "lastupdated", j[[1]]$lastupdated),
           file.path(EXT, "fetch_log.txt"))
