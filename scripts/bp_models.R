library(bpmodels)
library(tidyverse)
library(data.table)
library(doParallel)
library(foreach)
library(tictoc)

## Functions to reseed bp that are too small
min_tree_quick <- function(min_size = 100, stop_threshold = 100, 
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
  
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  
  chains <- foreach(
    chain=1:nrow(small_chains),
    .export=c("R","k","si","max_size","mlog","sdlog"),
    .packages=c("bpmodels", "tidyverse"),
    .combine="rbind"
  ) %dopar% {
    
    working_chain <- chains[chains$n == small_chains$n[chain], ]
    
    #only run until threshold is met
    while( 
      (working_chain %>% count(n) %>% pull(nn) %>% min()) < min_size
    ){
      if(chain_times == "non-overlap"){
        t_start <- max(working_chain$time) + 1
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
      working_chain <- rbind(working_chain, one_sim)
    }
    return(working_chain)
  }
  
  stopCluster(cl)
  
  return(chains)
}

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
    
    #this is where chains that have already met the threshold are
    # filtered out
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
# with small R, this can take a long (!) time to reach desired outbreak size (lot of processes die out)
n <- 10000
## R
R <- c(1.5) #, 2, 3)
## overdispersion
k <- c(0.35) # , 0.58, 1.18) # as per 0.58, 95% CI 0.35, 1.18). this is in the Bi et al medrXiv paper
## maximum outbreak size
max_size <- c(14962) #131071, 2391484)
## min_size we need the sum of branches in a scenario to be
min_size <- c(14962) #131071, 2391484)

Rk <- expand.grid(R = R, k = k)

# Serial interval
mlog <- log((4.7^2)/ (sqrt(2.9^2 + 4.7^2)))
vlog <- log(1 + (2.9^2/(4.7^2)))
sdlog <- sqrt(vlog) 

si <- function(x) round(rlnorm(x, meanlog = mlog, sdlog = sdlog))
trees <- list()

tictoc::tic()
for(i in 1:nrow(Rk)){
  Rk_i = Rk[i, ]
  Ri <- Rk_i$R
  ki <- Rk_i$k
  ## negative binomial offspring
  tree <- min_tree_quick(n=n, 
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
tictoc::toc()

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
  mutate(time_first_case_death = max(0, time - time_until_death)) %>% 
  select(-time)

tr_first_case_death
tr <- left_join(tr, tr_first_case_death)


tr_dt <- data.table::as.data.table(tr)
class(tr_dt$csum) <- "integer"
setorder(tr_dt, n, time)
tr_dt <- data.table::dcast(tr_dt, n + time_first_case_death ~ time, value.var="cid", fun.aggregate=sum)
tr_dt <- melt(tr_dt, id.vars=c("n","time_first_case_death"),variable.name="time",value.name="cases",variable.factor=FALSE, value.factor=FALSE)
class(tr_dt$time) <- "integer"
tr_dt[, cumcases := cumsum(cases), by="n"]

ggplot(
  tr_dt[n %in% sample(c(1:1000), 100, FALSE),],
  aes(x=time,y=cumcases, group=n)
)+geom_line(
  alpha=0.2
)+theme_bw(
)

## Difference in number of cases now and when the first death became a case

# when I ask for 1000 simulations, I (in this instance) get length(delta) == 1559
# I think we want this to be 1000 as well

# I think delta should be until simulated number of cases at current date minus delay
# onset-death
#delta <- max_case_death - tr$csum[tr$time == tr$time_first_case_death]
delta <- max_case_death - tr_dt[time == time_first_case_death, cumcases]

delta
hist(delta)

#P(delta)
table(delta)/length(delta)

delta_unique_p <- table(delta)/length(delta)
delta_p <- sapply(delta, function(x){delta_unique_p[as.character(x)]})

deaths <- 1
cfr <- 0.02
#P(deaths | delta)
#wonder whether to use cumulative number of deaths given cumulative number of
# cases at last case-observation minus delay of onset-to-death
dbinom(deaths, delta, cfr)

#P(deaths) ????
# P one death = cfr
# P x deaths = cfr^deaths ???
cfr^deaths

#P(delta | deaths)
# depending on k
dbinom(deaths, delta, cfr) * delta_p / cfr^deaths



## NUmber for the paper. 


# Max number paper
qgeom(0.975, prob = 0.005) # = 735
qgeom(0.975, prob = 0.01) # = 367


# max time until death

qgamma(0.975, shape = 4.726, rate = 0.3151) # 31.23776
hist(rgamma(10000, shape = 4.726, rate = 0.3151))

