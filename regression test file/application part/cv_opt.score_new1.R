# Arguments:
# Y:          n by K indicator matrix of classes (response)
# X:          n by s design matrix (covariates)
# R0:         s by s symmetric mass matrix
# R1:         s by s symmetric stiffness matrix
# kfolds:     the number of folds
# lambdas:    a set of tuning parameters of smoothness penalty
# ------------------------------------------------------------------------
# Outputs:
# beta:       s by K-1 matrix beta
# min_lambda: the tuning parameter which results in minimum error 
# cv_mse:     the error rate

cv_opt.score <- function(Y, X, R0, R1, kfolds, lambdas){
  
  K <- ncol(Y)
  
  nlambda <- length(lambdas)
  
  cv_mse <- rep(0, nlambda)
  
  fold <- sample(1:5, nrow(X), replace = TRUE)
  
  R0til.inv <- Matrix::Diagonal(x=1/Matrix::rowSums(R0))
  S<-R1 %*% sqrt(R0til.inv)
  for(i in 1:nlambda){
    print(i)
    error=0
    for(j in 1:kfolds){
      X_train <- X[!fold==j,]
      Y_train <- Y[!fold==j,]
      X_valid <- X[fold==j,]
      Y_valid <- Y[fold==j,]
      model <- opt.score(Y_train, X_train, R0, R1, lambdas[i])
      beta<-model$beta
      mean_fac<-model$mean_fac
      D <- 1/nrow(Y_valid) * (t(Y_valid)%*%Y_valid) # K by K matrix
      
      
      #L <- R1 %*% R0til.inv %*% R1
      #S <- Matrix::chol(L)
      #centered.X<-scale(X_valid,center=TRUE,scale=FALSE)
      
      centered.X<-sapply(1:length(mean_fac),function(x) X_valid[,x]-mean_fac[x])
      
      scoredX <- matrix(rep(0, nrow(Y_valid)*(K-1)), ncol=K-1)
      
      for(k in 1:(K-1)){
        scoredX[,k] <- as.vector(centered.X %*% beta[,k])        
      }
      
      ## LDA with two discriminant directions
      mu_k <- matrix(0, nrow = K-1, ncol = K) # centroid vectors
      pi_k <- vector(length = K) # priori probabilities
      sigma <- matrix(0, nrow = K-1, ncol = K-1) # covariance matrix
      
      for (l in 1:K){
        tmp <- matrix(scoredX[Y_valid[,l]==1,],ncol=K-1)
        pi_k[l] <- nrow(tmp)/nrow(Y_valid)
        mu_k[,l] <- colMeans(tmp)
        for(j in 1:nrow(tmp)){
          sigma <- sigma + (tmp[j,] - mu_k[,l]) %*% t(tmp[j,] - mu_k[,l])
        }
      }
      
      sigma <- sigma*(1/(nrow(Y_valid)-K))
      inv.sigma <- solve(sigma)
      
      deltatrain <- matrix(0, ncol = K, nrow = nrow(Y_valid))
      for (t in 1:K){
        deltatrain[,t] <- scoredX %*% inv.sigma %*% mu_k[,t]
        deltatrain[,t] <- deltatrain[,t] - 0.5*as.numeric(t(mu_k[,t]) %*% inv.sigma %*% mu_k[,t]) + log(pi_k[t])
      }
      
      yhattrain <- apply(deltatrain, 1, which.max)
      Y_valid_label <- apply(Y_valid, 1, which.max)
      error <- error + mean(yhattrain!=Y_valid_label) # training error rate?
    }
    cv_mse[i] <- error/kfolds
  }
  
  min_lambda <- as.numeric(lambdas[which(cv_mse==min(cv_mse))])
  beta <- opt.score(Y_train, X_train, R0, R1, min_lambda)
  
  return(list(beta=beta, min_lambda=min_lambda, cv_mse=cv_mse))
}
