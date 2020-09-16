# Arguments:
# Y:            n by 1 vector of the response variable.
# X:            n (number of rows) by p (number of columns) data matrix.
# R0:           p by p symmetric matrix of mass matrix.
# R1:           p by p symmetric matrix of stiff matrix.
# lambda:     tuning parameters of the penalty.
# ----------------------------------------------------------------------------
# Outputs:
# intercept:    intercept of the linear regression model.
# beta:         regression coefficients (slopes) of the linear regression model.
# ---------------------------------------------------------------------------

# lpl_penalty(Y=y,X=fdata,L=lpl,lambda.L=100,lambda.1=0.0001)
library(Matrix) 

lpl_penalty <- function(Y, X, R0,R1, lambda){


  ori.Y <- Y
  ori.X <- X
  
  # ----------------------
  # | Data Preprocessing |
  # ----------------------
  
  Y <- Y - mean(Y)
  n <- nrow(X)
  p <- ncol(X)
  scale.fac <- attr(scale(X), "scaled:scale")
  X <- scale(X)
  X[is.na(X)] <- 0
  
  # ----------------------
  # | create linear equation|
  # ----------------------
  # here we replace R1^-1 with a sparse diagonal matrix
  R1_inv=.sparseDiagona(1/rowsum(R1))
  L=R0*R1_inv*R0
  
  S<-chol(L)
  Xstar <- rbind(X, sqrt(lambda)*t(S)) 
  Ystar <- c(Y, rep(0, p))
  
  # ----------------------
  # | use glmnet to solve the euqation|
  # ----------------------
  
  betahat <- glmnet(Xstar, Ystar, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta
  beta <- betahat / scale.fac
  intercept <- mean(ori.Y - ori.X %*% truebetahat)
  return(list(intercept =   intercept, beta = beta))
}
