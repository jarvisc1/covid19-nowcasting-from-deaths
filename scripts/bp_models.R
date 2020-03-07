library(bpmodels)
library(tidyverse)

## Functions to reseed bp that are too small
min_tree <- function(min_size = 100, stop_threshold = 100, ..., 
                     chain_times = c("same", "non-overlap", "interval"), interval = 0) {
  ## generate chains
  chains <- chain_sim(infinite=stop_threshold, tree=T, ...)
  
  ## Find chains that are small
  small_chains <- chains %>%
    count(n) %>%
    filter(nn < min_size)
  
  # Check if all chains have taken off
  min_chain_n <- min(small_chains$nn)
  
  # Put start time for bp
  t_start = 0
  
  # Keep running bp untill all are at min size
  while(min_chain_n < min_size){
    
    small_chains <- chains %>%
      group_by(n) %>%
      summarise(nn = n(),
                t_max = max(time)) %>% 
      filter(nn < min_size)
    
   if(chain_times == "non-overlap"){
    t_start <- small_chains$t_max[chain] + 1
    print(t_start)
   } else if(chain_times == "interval"){
    t_start <- t_start + interval
  }
    
    
    for (chain in 1:nrow(small_chains)){
      one_sim <- chain_sim(1,infinite = max_size, 
                           offspring="nbinom", 
                           mu = R, 
                           size = k, 
                           tree = T, 
                           serial = si,
                           t0 = t_start)
      
      one_sim$n <- small_chains$n[chain]
      chains <- rbind(chains, one_sim)
    }
    min_chain_n <- chains %>%
      count(n) %>% 
      pull(nn) %>% 
      min()
  }
  return(chains)
}

## number of outbreaks to simulate
n <- 100
## R
R <- c(1.5, 2, 3)
## overdispersion
k <- c(0.35, 0.58, 1.18) # as per 0.58, 95% CI 0.35, 1.18). this is in the Bi et al medrXiv paper
## maximum outbreak size
max_size <- 2000
## min_size
min_size <- 2000

Rk <- expand.grid(R = R, k = k)

# Serial interval
mlog <- log((4.7^2)/ (sqrt(2.9^2 + 4.7^2)))
vlog <- log(1 + (2.9^2/(4.7^2)))
sdlog <- sqrt(vlog) 

si <- function(x) round(rlnorm(x, meanlog = mlog, sdlog = sdlog))


trees <- list()
for(i in 1:nrow(Rk)){
  Rk_i = Rk[i, ]
  Ri <- Rk_i$R
  ki <- Rk_i$k
  ## negative binomial offspring
  tree <- min_tree(n=n, 
                   min_size = min_size,
                   stop_threshold = max_size, 
                   offspring="nbinom", 
                   mu = R, 
                   size = k, 
                   serial = si)
  
  trees[[i]] <- list(Rk = Rk_i, bp_trees = tree)
  
}


# Max number paper
qgeom(0.975, prob = 0.04) # = 735
qgeom(0.975, prob = 0.01) # = 367


# max time until death

qgamma(0.975, shape = 4.726, rate = 0.3151) # 31.23776
hist(rgamma(10000, shape = 4.726, rate = 0.3151))

