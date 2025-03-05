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


# 新的 samples.iCAR 函数：生成 tau 和 log.risk 的一一对应列表
generate_tau_logrisk <- function(Data=NULL, tau.range=c(4, 400), intercept=0, k=100, s=100){

  # 定义 tau 范围，生成 k 个 tau 值
  tau <- seq(tau.range[1], tau.range[2], length.out=k)

  # 创建一个空的列表，用于存储结果
  result <- list()

  # 遍历每个 tau 值，生成对应的 log.risk
  for (i in 1:k) {
    # 获取当前的 tau 值
    current_tau <- tau[i]

    # 生成随机效应 xi，通过 rMVNormP_eigen
    xi <- rMVNormP_eigen(s, mu=rep(0, nrow(Rs)), Q=current_tau * Rs, Q.eigen=Rs.eigen)

    # 计算 log.risk
    log_risk <- intercept + xi  # 对于每个 tau，生成 s 个 log.risk 样本

    # 将当前 tau 和对应的 log.risk 组合成一行
    result[[i]] <- cbind(rep(current_tau, s), log_risk)  # 每行是 [tau, log_risk]
  }

  # 将所有结果合并为一个数据框（或者列表）
  result_df <- do.call(rbind, result)  # 合并所有 tau 和 log_risk
  colnames(result_df) <- c("tau", "log_risk")  # 设置列名为 tau 和 log_risk

  return(result_df)
}

