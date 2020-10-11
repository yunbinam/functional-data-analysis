# Arguments:
# Y:            n by 1 vector of the response variable.
# X:            n (number of rows) by p (number of columns) data matrix.
# R0:           p by p symmetric matrix of mass matrix.
# R1:           p by p symmetric matrix of stiff matrix.
# kfolds:       the folds number of cross validation.
# lambdas:      a sequence of tuning parameters of the smoothness penalty.
# ----------------------------------------------------------------------------
# Outputs:
# coef:         coeffcients of regression model, including intercept and beta.
# min_lambda:   the lambda which causes minimum cross validation MSE
# cv_mse:          the MSE of the best lambda on the whole dataset
# ---------------------------------------------------------------------------

##lambda is chosen by minimizing least squares of the regression (CV needed)
cv_smooth_ols_x<- function(Y, X, R0,R1, kfolds, lambdas){
  n<-nrow(X)
  p <- ncol(X)
  nlambda=length(lambdas)
  
  cv_mse=rep(0,nlambda)
  fold=cut(seq(1,nrow(X)),breaks=kfolds,labels=FALSE)
  
  ##dealing with matrix
  R0_inv<-Diagonal(x=1/Matrix::rowSums(R0))
  L<-R1%*%R0_inv%*%R1
  S<-chol(L)
  identity<-Diagonal(x=rep(1,p))
  
  
  for(i in 1:nlambda){
    error<-0
    lambda=lambdas[i]
    new_X<-matrix(0,nrow=n,ncol=p)
    Xstar <- rbind(identity, sqrt(lambda)*t(S))
    for(k in 1:nrow(X)){
      x_obs <- c(X[k,], rep(0, p))
      new_X[k,]<- as.vector(glmnet(Xstar, x_obs, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta)
    }
    
    
    for(j in 1:kfolds){
      
      X_train=new_X[!fold==j,]
      y_train=Y[!fold==j]
      X_valid=new_X[fold==j,]
      y_valid=Y[fold==j]
      
      # build functional regression model for smoothed X
      X_train=X_train%*%R0
      ori.y_train=y_train
      ori.X_train=X_train
      center_X_train <- scale(X_train,center=TRUE,scale=FALSE)
      center_y_train <- y_train - mean(y_train)
      
      
      beta=glmnet(center_X_train, center_y_train, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta
      intercept <- mean(ori.y_train - ori.X_train %*% beta)
      predict=intercept+X_valid%*%R0%*%beta
      error=error+mean((predict-y_valid)^2)
    }
    cv_mse[i]=error/kfolds
  }
  
  min_lambda=lambdas[which(cv_mse==min(cv_mse))]
  
  ## using selected lambda to build model based on all data----
  
  Xstar <- rbind(identity, sqrt(min_lambda)*t(S))
  for(k in 1:n){
    x_obs <- c(X[k,], rep(0, p))
    new_X[k,]<- as.vector(glmnet(Xstar, x_obs, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta)
  }
  
  X=new_X%*%R0
  ori.y=Y
  ori.X=X
  center_X <- scale(X,center=TRUE,scale=FALSE)
  center_y <- Y - mean(y)
  
  beta=glmnet(center_X, center_y, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta
  intercept <- mean(ori.y - ori.X %*% beta)
  coef=list(intercept=intercept,beta=beta)
  #coef=smooth_ols(Y, X, R0,R1, min_lambda)
  #intercept=coef$intercept
  #beta=coef$beta
  return(list(coef=coef,min_lambda=min_lambda,cv_mse=cv_mse))
}