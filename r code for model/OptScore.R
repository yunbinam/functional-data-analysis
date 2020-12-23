OptScore <- function(X_train,Y_train, R0, R1, lambda){
  
  Y <- as.matrix(Y_train)
  X <- as.matrix(X_train)
  n <- nrow(Y_train)
  K <- ncol(Y_train)
  s <- ncol(X_train)
  
  D <- 1/n * (t(Y_train)%*%Y_train) # K by K matrix
  R0til.inv <- Matrix::Diagonal(x=1/Matrix::rowSums(R0))
  S<-R1 %*% sqrt(R0til.inv)
  scale.fac <- attr(scale(X_train), "scaled:center")
  centered.X_train<-scale(X_train,center=TRUE,scale=FALSE) ##This centering is essential for the trivial solution to disappear
  Xstar <- rbind(centered.X_train, sqrt(lambda)*t(S))
  Q <- as.matrix(rep(1, K))
  beta <- matrix(rep(0, s*(K-1)), ncol=K-1)
  scoredX <- matrix(rep(0, n*(K-1)), ncol=K-1)
  
  for(i in 1:(K-1)){
    random.theta <- runif(K)
    theta <- (diag(K) - Q%*%t(Q)%*%D) %*% random.theta
    norm <- as.numeric(t(theta)%*%D%*%theta)
    norm.theta <- theta * sqrt(1/norm)
    mapped.Y <- Y_train %*% norm.theta
    Ystar <- c(mapped.Y, rep(0, s))
    beta[,i] <- as.numeric(glmnet::glmnet(Xstar, Ystar, 
                                          lambda=0, alpha=0, intercept=FALSE, 
                                          standardize=FALSE, thresh=1e-7)$beta)
    scoredX[,i] <- as.vector(centered.X_train %*% beta[,i])
    new.theta <- (diag(K) - Q%*%t(Q)%*%D) %*% solve(D)%*%t(Y)%*%X%*%beta[,i]
    new.norm <- as.numeric(t(as.matrix(new.theta))%*%D%*%as.matrix(new.theta))
    new.norm.theta <- new.theta * sqrt(1/new.norm)
    Q <- cbind(Q, new.norm.theta)
  }
  
  scoredX <- matrix(rep(0, n*(K-1)), ncol=K-1)
  
  for(k in 1:(K-1)){
    scoredX[,k] <- as.vector(centered.X_train %*% beta[,k])        
  }
  
  mu_k <- matrix(0, nrow = K-1, ncol = K) # centroid vectors
  pi_k <- vector(length = K) # priori probabilities
  sigma <- matrix(0, nrow = K-1, ncol = K-1) # covariance matrix
  
  for (l in 1:K){
    tmp <- matrix(scoredX[Y_train[,l]==1,],ncol=K-1)
    pi_k[l] <- nrow(tmp)/nrow(Y_train)
    mu_k[,l] <- colMeans(tmp)
    for(z in 1:nrow(tmp)){
      sigma <- sigma + (tmp[z,] - mu_k[,l]) %*% t(tmp[z,] - mu_k[,l])
    }
  }
  
  sigma <- sigma/(n-K)
  inv.sigma <- solve(sigma)
  deltatrain <- matrix(0, ncol = K, nrow = n)
  
  for (t in 1:K){
    deltatrain[,t] <- scoredX %*% inv.sigma %*% mu_k[,t]
    deltatrain[,t] <- deltatrain[,t] - 0.5*as.numeric(t(mu_k[,t]) %*% inv.sigma %*% mu_k[,t]) + log(pi_k[t])
  }
  
  model = structure(list(beta=beta,inv.sigma =inv.sigma , mu_k= mu_k,pi_k=pi_k,mean_fac=scale.fac), 
                    class = "OptScore")
  return(model)
}



predict.OptScore<-function(model,X_test,Y_test){
  
  mean_fac<-model$mean_fac
  centered.X_test<-sapply(1:length(mean_fac),function(x) X_test[,x]-mean_fac[x])
  mu_k <- model$mu_k # centroid vectors
  pi_k <- model$pi_k # priori probabilities
  inv.sigma <- model$inv.sigma # covariance matrix
  K=ncol(model$mu_k)
  beta=model$beta
  scoredX<-centered.X_test%*%beta
  deltatrain <- matrix(0, ncol = K, nrow = nrow(Y_test))
  
  for (t in 1:K){
    deltatrain[,t] <- scoredX %*% inv.sigma %*% mu_k[,t]
    deltatrain[,t] <- deltatrain[,t] - 0.5*as.numeric(t(mu_k[,t]) %*% inv.sigma %*% mu_k[,t]) + log(pi_k[t])
  }
  
  yhat <- apply(deltatrain, 1, which.max)
  Y_label <- apply(Y_test, 1, which.max)
  accuracy<- mean(yhat==Y_label)
  return(accuracy)
}




cv.OptScore <- function( X,Y, R0, R1, kfolds, lambdas){
  
  K <- ncol(Y)
  nlambda <- length(lambdas)
  cv_error <- rep(0, nlambda)
  fold <- sample(1:kfolds, nrow(X), replace = TRUE)
  R0til.inv <- Matrix::Diagonal(x=1/Matrix::rowSums(R0))
  S<-R1 %*% sqrt(R0til.inv)
  
  for(i in 1:nlambda){
    error=0
    for(j in 1:kfolds){
      X_train <- X[!fold==j,]
      Y_train <- Y[!fold==j,]
      X_valid <- X[fold==j,]
      Y_valid <- Y[fold==j,]
      model <- OptScore(X_train,Y_train, R0, R1, lambdas[i])
      err<-1-predict(model,X_valid,Y_valid)
      error <- error + err # training error rate?
    }
    cv_error[i] <- error/kfolds
  }
  
  min_lambda <- as.numeric(lambdas[which(cv_error==min(cv_error))])
  return(list(min_lambda=min_lambda, cv_error=cv_error))
}
