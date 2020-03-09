suppressPackageStartupMessages({
  require(jsonlite)
})

.args <- c("parameters.json", "parameters.rda")

#' the mean, sd of X ~ lognormal
#' does not directly correspond to the mean, sd of the related
#' normally distributed value, which is what (dpqr)lnorm uses for
#' parameters
lnorm_transform <- function(lnormmu, lnormsd) {
  return(list(
    meanlog = 2*log(lnormmu) - 0.5*log(lnormmu^2+ lnormsd^2),
    sdlog = sqrt(log(1+(lnormsd/lnormmu)^2))
  ))
}

nbinom_transform <- function(nbmu, nbdispersion) list(
  mu=nbmu, size=nbdispersion
)

.pars <- read_json(.args[1], simplifyVector = TRUE)

.rgen <- function(
  type, pars, f = switch(type,
    gamma = rgamma, lognormal = rlnorm,
    geom = rgeom, negbinomial = rnbinom
  )
) function(n) do.call(f, c(list(n=n), pars))

.qgen <- function(
  type, pars, f = switch(type,
    gamma = qgamma, lognormal = qlnorm,
    geom = qgeom, negbinomial = qnbinom)
) function(p) do.call(f, c(list(p=p), pars))

#' these functions are just draws, since parameters set them
rserial <- with(.pars$serial, .rgen(type, lnorm_transform(mu, sd)))
qserial <- with(.pars$serial, .qgen(type, lnorm_transform(mu, sd)))

rdeath <- with(.pars$death, .rgen(type, lnorm_transform(mu, sd)))
qdeath <- with(.pars$death, .qgen(type, lnorm_transform(mu, sd)))

#' this allows for different scenario inputs
roffspring <- with(.pars$offspring,
  function(n, R, k, pars = nbinom_transform(R, k),
    f=switch(type,
      gamma = rgamma, lognormal = rlnorm,
      geom = rgeom, negbinomial = rnbinom
  )) do.call(f, c(n=n, pars))
)

save(list=ls(), file=tail(.args, 1))