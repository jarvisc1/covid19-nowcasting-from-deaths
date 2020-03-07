library(bpmodels)
library(tidyverse)

## Functions to reseed bp that are too small
min_tree <- function(min_size = 100, stop_threshold = 100, 
                     chain_times = c("same", "non-overlap", "interval"), 
                     interval = 0, ...) {
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
    
    for (chain in 1:nrow(small_chains)){
      
      if(chain_times == "non-overlap"){
       t_start <- small_chains$t_max[chain] + 1
      } else if(chain_times == "interval"){
       t_start <- t_start + interval
      }
     
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


## Small values for testing
## number of outbreaks to simulate
n <- 10
## R
R <- c(1.5) #, 2, 3)
## overdispersion
k <- c(0.35) # , 0.58, 1.18) # as per 0.58, 95% CI 0.35, 1.18). this is in the Bi et al medrXiv paper
## maximum outbreak size
max_size <- 100
## min_size
min_size <- 100

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
                   serial = si,
                   chain_times = "non-overlap",
                   )
  
  trees[[i]] <- list(Rk = Rk_i, bp_trees = tree)
  
}


# Get one scenario
#####
tr <- trees[[1]]
tr <- tr$bp_trees 
tr %>% count(n)
tr <- tr %>% arrange(n, time)
 

# Count cumulative cases
tr$cid <- 1
tr$csum <- ave(tr$cid, tr$n, FUN = cumsum)

tr <- tr %>% group_by(n) %>% 
  mutate(csum_max = max(csum))

## Check non-overlap works
tr %>% filter(id == 1)


## Max cases for a death
max_case_death <- qgeom(0.975, prob = 0.005) # = 735
max_case_death <- 20 # test valye

# Tag in the data when this happens
tr <- tr %>% 
  mutate(above = if_else(csum >= max_case_death, 1, 0),
         time_step = if_else(above != lag(above) & above == 1, 1, 0)
  )


## Time until death 
time_until_death <- qgamma(0.975, shape = 4.726, rate = 0.3151) # 31.23776
time_until_death <- 3 #  Test for 2 days


tr_first_case_death <- tr[tr$time_step==1,] %>% 
  select(n, time) %>% 
  mutate(time_first_case_death = if_else(time - time_until_death < 0 ,
                                         0, 
                                         time - time_until_death)) %>% 
  select(-time)

tr_first_case_death
tr <- left_join(tr, tr_first_case_death)


## Difference in number of cases now and when the first death became a case
delta <- max_case_death - tr$csum[tr$time == tr$time_first_case_death]

delta


hist(delta)




## NUmber for the paper. 


# Max number paper
qgeom(0.975, prob = 0.005) # = 735
qgeom(0.975, prob = 0.01) # = 367


# max time until death

qgamma(0.975, shape = 4.726, rate = 0.3151) # 31.23776
hist(rgamma(10000, shape = 4.726, rate = 0.3151))

