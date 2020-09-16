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
# MSE:          the MSE of the best lambda on the whole dataset
# ---------------------------------------------------------------------------

library(Matrix) 

cv_spareg <- function(Y, X, R0,R1, kfolds, lambdas){
  
  nlambda=length(lambdas)
  
  cv_mse=seq(0,nlambda)
  fold=cut(seq(1,nrow(X)),breaks=kolds,labels=FALSE)
  
  for(i in 1:nlambda){
    error=0
    for(k in 1:kfolds){
      x_train=x[!fold==i,]
      y_train=y[!fold==i,]
      x_valid=x[fold==i,]
      y_valid=y[fold==i,]
      coef=spareg(y_train, x_train, R0,R1, lambdas[i])
      intercept=coef$intercept
      beta=coef$beta
      predict=intercept+x_valid%*%beta
      error=error+mean((predict-y_valid)^2)
    }
    cv_mse[i]=error/k
  }
  
  
  min_lambda=lambdas[which(CV_err==min(CV_err))]
  
  coef=spareg(Y, X, R0,R1, min_lambda)
  intercept=coef$intercept
  beta=coef$beta
  predict=intercept+X%*%beta
  mse=mse+mean((predict-Y)^2)
  
  return(list(coef=coef,min_lambda=min_lambda,MSE=mse))
}