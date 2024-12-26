library("AMISforInfectiousDiseases")

# Define simple "transmission" model where prevalence equals first parameter
transmission_model_identity <- function(seeds, parameters, n_tims=1) {
  return(matrix(parameters[,1], ncol=1))
}
# Number of locations
L <- 3
# Number of map samples
M <- 1000
prevalence_map <- matrix(NA, L, M)
# Produce samples for prevalence map with 3 locations given by B(2,1), B(1,1)=Uniform, B(1,2). 
for (l in 1:L) {
  prevalence_map[l,] <- rbeta(M, max(1,l-1), max(1,3-l))
}
# 2D exponential prior
rprior <- function(n) {
  params <- matrix(NA, n, 2)
  colnames(params) <- c("a","b")
  params[,1] <- rexp(n)
  params[,2] <- rexp(n)
  return(params)
}
dprior <- function(x, log=FALSE) {
  if (log) {
    return(sum(dexp(x, log=TRUE)))
  } else {
    return(prod(dexp(x)))
  }
}
prior <- list(rprior=rprior,dprior=dprior)
amis_params <- default_amis_params()
output <- amis(prevalence_map, transmission_model_identity, prior, amis_params, seed=1)

print(output)
summary(output)

original_par <- par()
par(cex.lab=1.5, cex.main=1.5, mar=c(5,4.5,4,2)+0.1)

par(mfrow=c(1,2))
plot_mixture_components(output, what = "uncertainty", cex=3)
plot_mixture_components(output, what = "density", nlevels=200)

par(mfrow=c(3,3))
plot(output, what = "a", type="hist", locations=1:L, breaks = 100)
plot(output, what = "b", type="hist", locations=1:L, breaks = 100)
plot(output, what = "prev", type="hist", locations=1:L, time=1, breaks = 100)

par(mfrow=c(1,3))
plot(output, what=c("a","b","prev"), type="CI", locations=1:L, 
     cex=3, lwd=3, measure_central="median")

calculate_summaries(output, what="prev", locations=1:L, alpha=0.05)

par(original_par)
