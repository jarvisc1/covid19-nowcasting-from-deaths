library(bpmodels)

## number of outbreaks to simulate
n <- 10000
## R
R <- c(1.5, 2, 3)
## overdispersion
k <- c(0.35, 0.58, 1.18) # as per 0.58, 95% CI 0.35, 1.18). this is in the Bi et al medrXiv paper
## maximum outbreak size
max_size <- 10000

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
  tree <- chain_sim(n=n, 
                     infinite = max_size, 
                     offspring="nbinom", 
                     mu = R, 
                     size = k, 
                     tree = T, 
                     serial = si)
  
  trees[[i]] <- list(Rk = Rk_i, bp_trees = tree)
  
}


# Max number paper
qgeom(0.975, prob = 0.005) # = 735
