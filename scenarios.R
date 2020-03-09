suppressPackageStartupMessages({
  require(jsonlite)
  require(data.table)
})

.args <- c("scenarios.json", "scenarios.rds")
.args <- commandArgs(trailingOnly = TRUE)

pars <- read_json(.args[1], simplifyVector = T)

res <- data.table(do.call(expand.grid, pars))

saveRDS(res, tail(.args, 1))