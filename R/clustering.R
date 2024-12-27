#' Fit an MVN mixture model using mclust
#' 
#' Can handle 1D as well as multivariate clustering. \code{dat} must have \code{nrow} observation and \code{ncol}
#'  dimensions, even if the number of dimensions is 1.
#' Uses BIC to determine the best number of components, up to max.components.
#'
#' @param dat A MATRIX or dataframe containing the observations to cluster.
#' @param max.components A positive integer specifying the maximum number of components to fit in the mixture.
#' @return list containing selected output from \code{mclust}:
#' \describe{
#' \item{\code{G}}{the best number of components G.}
#' \item{\code{probs}}{vector of cluster probabilities (mixing weights).}
#' \item{\code{Mean}}{matrix of cluster means.}
#' \item{\code{Sigma}}{array of cluster covariance matrices.}
#' \item{\code{BIC}}{The BIC of the chosen mixture.}
#' \item{\code{ModelName}}{The model name from the package \code{mclust}.}
#' }
#' @import mclust
#' @noRd
# #' @export
fit_mixture <- function(dat,max.components=10) {
  n <- nrow(dat)
  d <- ncol(dat)
  colnames(dat) <- NULL # remove colnames to prevent instigating bug in mclust
  if (n<d+1) {stop("Not enough observations to fit mixture model.")}
  max.components <- min(max.components,floor(n/(d+1)))
  BICmatrix <- matrix(NA,nrow=max.components,ncol=1)
  # Start by fitting one group
  G <- 1 # number of groups
  if (d==1) {
    modelName <- "X" 
  } else {
    modelName <- "XXX"
  }
  clustering <- mclust::mvn(modelName=modelName,data=dat)
  BIC <- mclust::bic(modelName=modelName,loglik=clustering$loglik,n=n,d=d,G=1)
  message(paste0("  BIC: ", BIC))
  BICmatrix[1, 1] <- BIC
  # fit agglomerative clustering model
  if (d==1) {
    modelName <- "V"
  } else {
    modelName <- "VVV"
  }
  hcPairs <- mclust::hc(modelName=modelName,data=dat)
  cut.tree <- mclust::hclass(hcPairs=hcPairs,G=2:max.components)
  for (g in 2:max.components) {
    z <- mclust::unmap(cut.tree[,g-1]) # extract cluster indices
    # Run EM algorithm
    em <- mclust::me(modelName=modelName,data=dat,z=z)
    em$BIC <- mclust::bic(modelName=modelName,loglik=em$loglik,n=n,d=d,G=g)
    BICmatrix[g, 1] <- em$BIC
    if (!is.na(em$BIC) && em$BIC>BIC) { # mclust:::bic calculates the negative BIC
      clustering <- em
      G <- g
      BIC <- em$BIC
      clustering$z <- em$z
    }
  }
  colnames(BICmatrix) <- modelName
  rownames(BICmatrix) <- 1:max.components
  clustering$BIC <- BICmatrix
  clustering$data <- dat
  clustering$hypvol <- NA
  class(clustering) <- "Mclust"
  return(list(G=G,probs=clustering$parameters$pro,Mean=matrix(clustering$parameters$mean,d,G),
              Sigma=array(clustering$parameters$variance$sigma,c(d,d,G)),BIC=BIC,modelName=clustering$modelName, 
              clustering=clustering))
}
