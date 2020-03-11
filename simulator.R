suppressPackageStartupMessages({
  require(bpmodels)
  require(data.table)
})

consolidator <- function(dt, infinite, offset = 0) if (nrow(dt)) dt[,
  .SD[order(time)][1:min(.N,infinite)], by=n
][,
  .(cases = .N), keyby = .(n, time)
][,
  cumcases := cumsum(cases) + offset, by=n
] else data.table(n=integer(0), time=integer(0), cases=integer(0), cumcases=integer(0))

min_tree_quick <- function(
  target_cases, 
  chain_times = c("same", "non-overlap", "interval"), 
  interval = 0, debug = 0, tree = TRUE,
  rngseed = as.integer(Sys.time()), ...
) {
  set.seed(rngseed)
  if (debug > 1) browser()
  ## generate chains
  chains <- data.table(chain_sim(infinite=target_cases*1.5, tree=tree, ...))
  # don't need these; if behavior of bpmodels changes, will need to remove
  chains$id <- chains$ancestor <- chains$generation <- NULL
  
  small_chains <- chains[, if(.N < target_cases) .SD, by=n]
  big_chains <- consolidator(chains[,
    if (.N >= target_cases) .SD,
    by=n
  ], target_cases)
  
  # expected number of runs to get one take off
  # always do at least 2
  extra_runs <- ceiling(
    big_chains[,max(2, chains[,length(unique(n))]/length(unique(n))) ]
  )
  
  pars <- list(...)
  pars$n <- extra_runs
  pars$tree <- tree
  
  expanded_chains <- small_chains[,{
    subdt <- copy(.SD)[, n := .BY ]
    consolidated <- consolidator(subdt, target_cases)
    if (debug > 0) browser()
    while(consolidated[.N, cumcases <= target_cases ]) {
      reft <- consolidated[.N, time+1]
      offst <- consolidated[.N, cumcases]
      more_chains <- data.table(do.call(
        chain_sim, c(infinite=target_cases*1.5, pars)
      ))
      more_chains$id <- more_chains$ancestor <- more_chains$generation <- NULL
      mxtms <- cumsum(c(reft, head(more_chains[, max(time)+1, by = n]$V1, -1)))
      more_chains[,time := time + mxtms[.GRP], by=n]
      more_chains[, n := .BY ]
      consolidated <- rbind(
        consolidated,
        consolidator(more_chains, target_cases, offst)
      )
    }
    consolidated
  }, by=n]
  
  if (debug > 2) browser()
  
  return(rbind(big_chains, expanded_chains[, -1]))
}

#' @examples 
#' n <- 20
#' Ri <- 1.5
#' ki <- 0.35
#' need_size <- 14962L
#' si <- function(x,
#'   mlog=log((4.7^2)/ (sqrt(2.9^2 + 4.7^2))),
#'   sdlog=sqrt(log(1 + (2.9^2/(4.7^2))))
#' ) round(rlnorm(x, meanlog = mlog, sdlog = sdlog))
#' tree <- min_tree_quick(
#'   n=n, target_cases = need_size,
#'   offspring="nbinom", mu = Ri, size = ki,
#'   serial = si, chain_times = "non-overlap"
#'   , debug = 3, rngseed = 1234L
#' )
#' 
#' require(ggplot2)
#' p <- ggplot(tree) + aes(time, cumcases, group=n) +
#' geom_step(alpha = 0.1) + theme_minimal() +
#' geom_hline(yintercept = need_size, color = "red") +
#' scale_x_continuous("day") +
#' scale_y_continuous("cumulative cases")

#' 
#' 
#' 
#' 