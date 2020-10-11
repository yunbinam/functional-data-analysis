# Arguments:
# Y:            n by 1 vector of the response variable.
# X:            n (number of rows) by p (number of columns) data matrix.
# R0:           p by p symmetric matrix of mass matrix.
# R1:           p by p symmetric matrix of stiff matrix.
# lambdas:      a sequence of tuning parameters of the smoothness penalty.
# ----------------------------------------------------------------------------
# Outputs:
# coef:         coeffcients of regression model, including intercept and beta.
# min_lambda:   the lambda which causes minimum re-construction error
# ---------------------------------------------------------------------------

## lambda is chosen by minimizing the reconstruction error of X (no cv neeeded)
smooth_recon_x<- function(Y, X, R0,R1, lambdas){
  n<-nrow(X)
  p <- ncol(X)
  nlambda=length(lambdas)
  
  mse=rep(0,nlambda)
  
  ##dealing with matrix
  R0_inv<-Diagonal(x=1/Matrix::rowSums(R0))
  L<-R1%*%R0_inv%*%R1
  S<-chol(L)
  identity<-Diagonal(x=rep(1,p))
  
  stored_X<-matrix(0,nrow=n,ncol=p)
  new_X<-matrix(0,nrow=n,ncol=p)
  min_error<-10000
  for(i in 1:nlambda){
    error<-0
    lambda=lambdas[i]
    
    Xstar <- rbind(identity, sqrt(lambda)*t(S))
    
    
    for(k in 1:n){
      x_obs <- c(X[k,], rep(0, p))
      new_X[k,]<- as.vector(glmnet(Xstar, x_obs, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta)
      error=error+mean((X[k,]-new_X[k,])^2)
    }
    if(min_error>error) 
    {
      stored_X<-new_X
      min_error=error
    }
    mse[i]=error
  }
  
  min_lambda=lambdas[which(mse==min(mse))]
  Xstar <- rbind(identity, sqrt(min_lambda)*t(S))
  new_X<-stored_X
  
  # for(k in 1:nrow(X)){
  #   x_obs <- c(X[k,], rep(0, p))
  #   new_X[k,]<- as.vector(glmnet(Xstar, x_obs, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta)
  # }
  
  
  ## using selected lambda to build model based on all data
  
  X=new_X%*%R0
  ori.y=Y
  ori.X=X
  center_X <- scale(X,center=TRUE,scale=FALSE)
  center_y <- Y - mean(y)
  
  beta=glmnet(center_X, center_y, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta
  intercept <- mean(ori.y - ori.X %*% beta)
  coef=list(intercept=intercept,beta=beta)
  
  return(list(coef=coef,min_lambda=min_lambda))
}