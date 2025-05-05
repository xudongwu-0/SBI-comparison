###########################################################################################################
## Auxiliary function to sample from a singular multivariate Normal using its precision parameterization ##
###########################################################################################################
rMVNormP_eigen <- function(s, mu, Q, Q.eigen=NULL){
  
  p <- length(mu)
  
  # Compute eigen decomposition
  if(is.null(Q.eigen)) Q.eigen <- eigen(Q,symmetric=TRUE)
  
  # Small negative eigenvalues can occur just due to numerical error and should
  # be set back to zero. 
  Q.eigen$values[(Q.eigen$values < 1e-12) & (Q.eigen$values > -1e-12)] <- 0
  
  # If any eigen values are large negative throw error (something wrong)
  if (any(Q.eigen$values < -1e-12)) stop("Non-trivial negative eigenvalues present")
  
  # Calculate square root and reveal rank (k)
  k <- sum(Q.eigen$values > 0)
  pos <- which(Q.eigen$values > 0)
  L <- Q.eigen$vectors[,pos] %*% Diagonal(k,1/sqrt(Q.eigen$values[pos]))
  
  # Sample from Z and transform them into samples of X 
  Z <- as(matrix(rnorm(k*s), k, s),"Matrix")
  X <- L%*%Z
  X <- as(sweep(X, 1, mu, FUN=`+`),"matrix")
}


###########################################################################################
## Algorithm to generate training data for the disease mapping neural network            ##
##                                                                                       ##
## + Step 1: Generate "k" values of the spatial precision parameter on a regular grid    ##
## + Step 2: For each value, simulate "s" random vectors from an iCAR prior distribution ##
## + Step 3: Compute the linear predictor (log-risk) of the model                        ##
## + Step 4: Generate the vector of observed cases from a Poisson distribution           ##
###########################################################################################
samples.iCAR <- function(Data=NULL, tau.range=c(4,400), intercept=0, k=100, s=100, l=100){

  tau <- seq(tau.range[1],tau.range[2],length.out=k)

  Obs <- lapply(tau, function(tau){

    #cat(sprintf("tau: %d \n",tau))
    set.seed(123)

    xi <- rMVNormP_eigen(s, mu=rep(0,nrow(Rs)), Q=tau*Rs, Q.eigen=Rs.eigen)

    log.risk <- intercept + xi

    lapply(seq(l), function(l){
      apply(log.risk, 2, function(x) rpois(nrow(Rs), Data$E*exp(x)))
    })

    Observations <- lapply(seq(l), function(i){
      apply(log.risk, 2, function(x) rpois(nrow(Rs), Data$E * exp(x)))
    })

    # 返回包含 Observed_Cases、Log_Risk 和 Tau 的列表
    list(Observed_Cases = Observations)

  })
  names(Obs) <- paste0("tau=",tau)

  return(Obs)
}
