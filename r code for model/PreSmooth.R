PreSmooth<- function(X_train,Y_train,R0, R1,lambda){
  n<-nrow(X_train)
  p <- ncol(X_train)
  identity <- Diagonal(x=rep(1,p))
  R0til.inv <- Matrix::Diagonal(x=1/Matrix::rowSums(R0))
  S <- R1 %*% sqrt(R0til.inv)
  Xstar <- rbind(identity, sqrt(lambda)*t(S))
  smooth_X<- matrix(0,nrow=n,ncol=p)
  
  for(i in 1:n){
    x_obs <- c(X_train[i,], rep(0, p))
    smooth_X[i,] <- as.vector(glmnet(Xstar, x_obs, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta)
  }
  
  df.train <- data.frame(cbind(smooth_X, Y=Y_train))
  fit <- lda(Y~., df.train)
  fitted.y <- predict(fit, df.train)$class
  accuracy <- mean(Y_train == fitted.y)
  model = structure(list(fit=fit,lambda=lambda,accuracy=accuracy,R0=R0,R1=R1), 
                    class = "PreSmooth")
  return(model)
}



predict.PreSmooth<-function(model,X_test,Y_test,lambda){
  R0<-model$R0
  R1<-model$R1
  n<-nrow(X_test)
  p<-ncol(X_test)
  identity <- Diagonal(x=rep(1,p))
  R0til.inv <- Matrix::Diagonal(x=1/Matrix::rowSums(R0))
  S <- R1 %*% sqrt(R0til.inv)
  Xstar <- rbind(identity, sqrt(lambda)*t(S))
  smooth_X<- matrix(0,nrow=n,ncol=p)
  
  for(i in 1:n){
    x_obs <- c(X_test[i,], rep(0, p))
    smooth_X[i,] <- as.vector(glmnet(Xstar, x_obs, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta)
  }
  
  df.test <- data.frame(cbind(smooth_X, Y=Y_test))
  fit<-model$fit
  fitted.y <- predict(fit, df.test)$class
  accuracy <- mean(Y_test == fitted.y)
  return(accuracy)
}



cv.PreSmooth <- function(X,Y,R0, R1, folds,lambdas){
  
  nlambda<-length(lambdas)
  n<-nrow(X)
  p <- ncol(X)
  
  fold <- sample(1:folds, nrow(X), replace = TRUE)
  
  identity <- Diagonal(x=rep(1,p))
  R0til.inv <- Matrix::Diagonal(x=1/Matrix::rowSums(R0))
  S <- R1 %*% sqrt(R0til.inv)
  
  cv_error <- rep(0, nlambda)
  
  for(i in 1:nlambda){
    lambda=lambdas[i]
    Xstar <- rbind(identity, sqrt(lambda)*t(S))
    smooth_X<- matrix(0,nrow=n,ncol=p)
    
    err<-0
    for(j in 1:n){
      x_obs <- c(X[j,], rep(0, p))
      smooth_X[j,] <- as.vector(glmnet(Xstar, x_obs, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta)
    }
      for(k in 1:folds){
        X_train <- smooth_X[!fold==k,]
        Y_train <- Y[!fold==k]
        X_valid <- smooth_X[fold==k,]
        Y_valid <- Y[fold==k]
        
        df.train <- data.frame(cbind(X_train, Y=Y_train))
        df.valid <- data.frame(cbind(X_valid, Y=Y_valid))
        fit <- lda(Y~., df.train)
        fitted.y <- predict(fit, df.valid)$class
        error <- mean(Y_valid != fitted.y)
        
        err<-err+error
      }
    cv_error[i]<-err/folds
  }
  
  min_lam<-lambdas[which.min(cv_error)]
  
  return(list(min_lambda=min_lam,cv_error=cv_error))
  
}

