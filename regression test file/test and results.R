# get accuracy
test_acc<-function(model,x_test,y_test){
  K=ncol(y_test)
  beta=model$beta
  mu_k <- model$mu_k # centroid vectors
  pi_k <- model$pi_k # priori probabilities
  inv.sigma <- model$inv.sigma # covariance matrix
  
  center.test<-sapply(1:length(model$mean_fac),function(x) x_test[,x]-model$mean_fac[x])
  scoredX<-center.test%*%beta
  
  deltatrain <- matrix(0, ncol = K, nrow = nrow(y_test))
  for (t in 1:K){
    deltatrain[,t] <- scoredX %*% inv.sigma %*% mu_k[,t]
    deltatrain[,t] <- deltatrain[,t] - 0.5*as.numeric(t(mu_k[,t]) %*% inv.sigma %*% mu_k[,t]) + log(pi_k[t])
  }
  
  yhat <- apply(deltatrain, 1, which.max)
  Y_label <- apply(y_test, 1, which.max)
  accuracy<- mean(yhat==Y_label)
  list(accuracy=accuracy,yhat=yhat)
}
#----------------------------
# example
# given x_train,y_train, x_test, y_test where y is indicator matrices
# 
# cv_model<-cv_opt.score(y_train, x_train, R0_642, R1_642, 5,10^seq(-5,5,1))
# model1<-opt.score(y_train, x_train, R0_642, R1_642, cv_model$min_lambda)
# result1<-test_acc(model1,x_test,y_test)
# result1$accuracy
# result1$yhat