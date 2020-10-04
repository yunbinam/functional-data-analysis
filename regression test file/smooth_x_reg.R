# Arguments:
# Y:            n by 1 vector of the response variable.
# X:            n (number of rows) by p (number of columns) data matrix.
# R0:            p by p symmetric matrix of the penalty weighted laplacian matrix.
# R1:     tuning parameters of the penalty weighted laplacian matrix.
# lambda:     tuning parameters of the smoothness penalty
# ----------------------------------------------------------------------------
# Outputs:
# intercept:    intercept of the linear regression model.
# beta:         regression coefficients (slopes) of the linear regression model.
# ---------------------------------------------------------------------------

smooth_x_reg <- function(Y,X, R0,R1, lambda){
  start.time <- Sys.time()
  n <- nrow(X)
  p <- ncol(X)
  new_X=matrix(0,nrow=n,ncol=p)
  identity<-Diagonal(x=rep(1,p))
  
  R0_inv<-Diagonal(x=1/Matrix::rowSums(R0))
  L<-R1%*%R0_inv%*%R1
  S<-chol(L)
  Xstar <- rbind(identity, sqrt(lambda)*t(S))
  for(i in 1:n){
    print(i)
    x_obs <- c(X[i,], rep(0, p))
    start.time <- Sys.time()
    new_X[i,]<- as.vector(glmnet(Xstar, x_obs, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta)
    end.time <- Sys.time()
    print(end.time-start.time)
  }
  
  input_X=new_X%*%R0
  
  # build functional regression model for smoothed X
  ori.Y<-Y
  ori.X<- input_X
  Y <- Y - mean(Y)
  input_X<-scale(input_X,center=TRUE,scale=FALSE)
  beta=glmnet(input_X, Y, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta
  intercept <- mean(ori.Y - ori.X %*% beta)
  return(list(intercept =intercept, beta = beta))
}  

