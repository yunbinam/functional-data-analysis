FPCA_LDA<-function(X_train,Y_train,FEMbasis,lambda,d){
  
  scale.fac <- attr(scale(X_train), "scaled:center")
  center.X_train<- scale(X_train,center=TRUE,scale=FALSE)
  FPCA_solution <- FPCA.FEM(datamatrix = center.X_train,FEMbasis = FEMbasis, lambda = lambda, nPC = d)
  scores<-FPCA_solution$scores
  df.train <- data.frame(cbind(scores, y=Y_train))
  fit <- lda(y~., df.train)
  fitted.y <- predict(fit, df.train)$class
  accuracy <- mean(Y_train == fitted.y)
  model = structure(list(fit=fit,lambda=lambda,FEMbasis=FEMbasis,mean_fac=scale.fac,accuracy=accuracy,components=FPCA_solution$loadings.FEM$coeff), 
                    class = "FPCA_LDA")
  return(model)
}



predict.FPCA_LDA<-function(model,X_test,Y_test){
  w<-model$components
  fit<-model$fit
  lambda=model$lambda
  FEMbasis=model$FEMbasis
  mean_fac<-model$mean_fac
  centered.X_test<-sapply(1:length(mean_fac),function(x) X_test[,x]-mean_fac[x])
  scores<-centered.X_test%*%w
  df.test <- data.frame(cbind(scores, y=Y_test))
  fitted.y <- predict(fit, df.test)$class
  accuracy <- mean(Y_test == fitted.y)
  return(accuracy)
}


cv.FPCA_LDA<-function(X,Y,FEMbasis,lambdas,kfolds,d){
  nlambda <- length(lambdas)
  cv_error <- rep(0, nlambda)
  fold <- sample(1:kfolds, nrow(X), replace = TRUE)
  for(i in 1:nlambda){
    err<-0
    lambda<-lambdas[i]
    for(j in 1:kfolds){
      X_train <- X[!fold==j,]
      Y_train <- Y[!fold==j]
      X_valid <- X[fold==j,]
      Y_valid <- Y[fold==j]
      model<-FPCA_LDA(X_train,Y_train,FEMbasis,lambda,d)
      error<-1-predict(model,X_valid,Y_valid)
      err<-err+error
    }
    cv_error[i]<-err/kfolds
  }
  min_lambda<-lambdas[which.min(cv_error)]
  return(list(min_lambda=min_lambda,cv_error=cv_error))
}


