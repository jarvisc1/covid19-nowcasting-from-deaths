suppressPackageStartupMessages({
  require(bpmodels)
  require(data.table)
})

consolidator <- function(dt, offset = 0) if (nrow(dt)) dt[,
  .(cases = .N), keyby = .(n, time)
][,
  cumcases := cumsum(cases) + offset, by=n
] else data.table(n=integer(0), time=integer(0), cases=integer(0), cumcases=integer(0))

min_tree_quick <- function(
  min_size = 100, stop_threshold = 100, 
  chain_times = c("same", "non-overlap", "interval"), 
  interval = 0, ...
) {
  
  ## generate chains
  chains <- data.table(
    chain_sim(infinite=stop_threshold, tree=T, ...)
  )
  
  small_chains <- chains[, if(.N < min_size) .SD, by=n]
  big_chains <- consolidator(chains[,
    if (.N >= min_size) .SD,
    by=n
  ])
  
  # expected number of runs to get one take off
  # always do at least 2
  extra_runs <- ceiling(
    big_chains[,max(2, chains[,length(unique(n))]/length(unique(n))) ]
  )
  
  pars <- list(...)
  pars$n <- extra_runs
  pars$infinite <- stop_threshold
  pars$tree <- T
  
  expanded_chains <- small_chains[,{
    subdt <- copy(.SD)[, n := .BY ]
    consolidated <- consolidator(subdt)
    while(consolidated[.N, cumcases < stop_threshold ]) {
      reft <- consolidated[.N, time+1]
      more_chains <- data.table(do.call(chain_sim, pars))
      mxtms <- cumsum(c(reft, head(more_chains[, max(time)+1, by = n]$V1, -1)))
      more_chains[,time := time + mxtms[.GRP], by=n]
      more_chains[, n := .BY ]
      consolidated <- rbind(consolidated, consolidator(more_chains, consolidated[.N, cumcases]))
    }
    consolidated
  }, by=n]
  
  return(rbind(big_chains, expanded_chains[, -1]))
}

#' @examples 
#' n <- 20
#' Ri <- c(1.5) 
#' ki <- c(0.35)
#' max_size <- c(14962)
#' min_size <- max_size
#' si <- function(x,
#'   mlog=log((4.7^2)/ (sqrt(2.9^2 + 4.7^2))),
#'   sdlog=sqrt(log(1 + (2.9^2/(4.7^2))))
#' ) round(rlnorm(x, meanlog = mlog, sdlog = sdlog))
#' tree <- min_tree_quick(
#'   n=n, min_size = min_size, stop_threshold = max_size,
#'   offspring="nbinom", mu = Ri, size = ki,
#'   serial = si, chain_times = "non-overlap"
#' )
#' 
#' require(ggplot2)
#' ggplot(tree) + aes(x=time, y=cumcases, color=factor(n)) + geom_step()
#' 
#' 
#' 
#' 