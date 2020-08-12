# Arguments:
# Y:            n by 1 vector of the response variable.
# X:            n (number of rows) by p (number of columns) design matrix.
# L:            p by p symmetric matrix of the penalty weighted laplacian matrix.
# lambda.L:     tuning parameters of the penalty weighted laplacian matrix.
# lambda.1:     tuning parameters of the L_1 penalty.
# ----------------------------------------------------------------------------
# Outputs:
# intercept:    intercept of the linear regression model.
# beta:         regression coefficients (slopes) of the linear regression model.
# ---------------------------------------------------------------------------

# lpl_penalty(Y=y,X=fdata,L=lpl,lambda.L=100,lambda.1=0.0001)

lpl_penalty <- function(Y, X, L, lambda.L, lambda.1 = 0){
  
  ori.Y <- Y
  ori.X <- X
  
  # ----------------------
  # | Data Preprocessing |
  # ----------------------
  
  Y <- Y - mean(Y)
  n <- nrow(X)
  p <- ncol(X)
  #scale.fac <- attr(scale(X), "scaled:scale")
  #X <- scale(X)
  #no scale because of na
  
  # --------------------------------------------------------
  #  see Li and Li (2008) 'Network-constrained regularization and variable selection 
  #                       for analysis of genomic data' for reference
  # --------------------------------------------------------
  
  S<-chol(L)
  Xstar <- rbind(X, sqrt(lambda.L)*t(S)) / sqrt(1+lambda.L)
  Ystar <- c(Y, rep(0, p))
  gammastar <- lambda.1 / sqrt(1+lambda.L) 
  
  betahatstar <- glmnet(Xstar, Ystar, lambda = gammastar, alpha=1,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta
  betahat <- betahatstar / sqrt(1+lambda.L)
  truebetahat <- betahat
  truealphahat <- mean(ori.Y - ori.X %*% truebetahat)
  return(list(intercept = truealphahat, beta = truebetahat))
}