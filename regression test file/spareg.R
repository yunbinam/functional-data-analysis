# Arguments:
# Y:            n by 1 vector of the response variable.
# X:            n (number of rows) by p (number of columns) data matrix.
# R0:            p by p symmetric matrix of the penalty weighted laplacian matrix.
# R1:     tuning parameters of the penalty weighted laplacian matrix.
# lambda.1:     tuning parameters of the smoothness penalty
# ----------------------------------------------------------------------------
# Outputs:
# intercept:    intercept of the linear regression model.
# beta:         regression coefficients (slopes) of the linear regression model.
# ---------------------------------------------------------------------------


spareg <- function(Y, X, R0,R1, lambda){
  
  
  ori.Y <- Y
  ori.X <- X
  
  # ----------------------
  # | Data Preprocessing |
  # ----------------------
  
  Y <- Y - mean(Y)
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X,center=TRUE,scale=FALSE)
  
  
  # ----------------------
  # | linear equation|
  # ----------------------
  # here we replace R0^-1 with a sparse diagonal matrix
  R0_inv<-Diagonal(x=1/Matrix::rowSums(R0))
  L<-R1%*%R0_inv%*%R1
  
  S<-chol(L)
  
  Xstar <- rbind(X, sqrt(lambda)*t(S)) 
  Ystar <- c(Y, rep(0, p))
  
  # ----------------------
  # | use glmnet to solve the equation|
  # ----------------------
  
  beta <- glmnet(Xstar, Ystar, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta
  intercept <- mean(ori.Y - ori.X %*% beta)
  return(list(intercept =intercept, beta = beta))
}