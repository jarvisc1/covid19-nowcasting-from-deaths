suppressPackageStartupMessages({
  require(jsonlite)
  require(data.table)
})

.args <- c("scenarios.json", "scenarios.ssv")
.args <- commandArgs(trailingOnly = TRUE)

pars <- read_json(.args[1], simplifyVector = T)

# for branching process runs, only want R, k grid
pars$CFR <- NULL

fwrite(do.call(expand.grid, pars), tail(.args, 1), sep = " ", col.names = FALSE)