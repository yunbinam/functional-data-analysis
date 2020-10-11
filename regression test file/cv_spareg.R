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


cv_spareg <- function(y, X, R0,R1, kfolds, lambdas){
  
  nlambda=length(lambdas)
  
  cv_mse=rep(0,nlambda)
  fold=cut(seq(1,nrow(X)),breaks=kfolds,labels=FALSE)
  
  for(i in 1:nlambda){
    error=0
    for(j in 1:kfolds){
      X_train=X[!fold==j,]
      y_train=y[!fold==j]
      X_valid=X[fold==j,]
      y_valid=y[fold==j]
      coef=spareg(y_train, X_train, R0,R1, lambdas[i])
      intercept=coef$intercept
      beta=coef$beta
      predict=intercept+X_valid%*%beta
      error=error+mean((predict-y_valid)^2)
    }
    cv_mse[i]=error/kfolds
  }
  
  
  min_lambda=lambdas[which(cv_mse==min(cv_mse))]
  coef=spareg(y, X, R0,R1, min_lambda)
  intercept=coef$intercept
  beta=coef$beta
  
  
  return(list(coef=coef,min_lambda=min_lambda,cv_mse=cv_mse))
}