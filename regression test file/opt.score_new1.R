# Arguments:
# Y:          n by K indicator matrix of classes (response)
# X:          n by s design matrix (covariates)
# R0:         s by s symmetric mass matrix
# R1:         s by s symmetric stiffness matrix
# lambda:     tuning parameter of smoothness penalty
# ------------------------------------------------------------------------
# Outputs:
# beta:       s by K-1 matrix beta

## simple example

opt.score <- function(Y, X, R0, R1, lambda){
  
  
  
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  
  n <- nrow(Y)
  K <- ncol(Y)
  s <- ncol(X)
  
  ## 1. prepare matrices
  D <- 1/n * (t(Y)%*%Y) # K by K matrix
  
  R0til.inv <- Matrix::Diagonal(x=1/Matrix::rowSums(R0))
  #L <- R1 %*% R0til.inv %*% R1
  #S <- Matrix::chol(L)
  S<-R1 %*% sqrt(R0til.inv)
  
  
  scale.fac <- attr(scale(X), "scaled:center")
  centered.X<-scale(X,center=TRUE,scale=FALSE) ##This centering is essential for the trivial solution to disappear
  X.star <- rbind(centered.X, sqrt(lambda)*t(S))
  
  ## 3. pre-specify theta1
  Q <- as.matrix(rep(0, ncol(Y)))
  beta <- matrix(rep(0, s*(K-1)), ncol=K-1)
  scoredX <- matrix(rep(0, n*(K-1)), ncol=K-1)
  
  ## 4. loop
  for(i in 1:(K-1)){
    ## a. compute theta
    random.theta <- runif(K)
    theta <- (diag(K) - as.matrix(Q[,i])%*%t(as.matrix(Q[,i]))%*%D) %*% random.theta
    norm <- as.numeric(t(theta)%*%D%*%theta)
    norm.theta <- theta * sqrt(1/norm)
    
    ## b. solve beta
    mapped.Y <- Y %*% norm.theta
    #centered.mapped.Y <- mapped.Y - mean(mapped.Y)
    Y.star <- c(mapped.Y, rep(0, s))
    beta[,i] <- as.numeric(glmnet::glmnet(X.star, Y.star, 
                                          lambda=0, alpha=0, intercept=FALSE, 
                                          standardize=FALSE, thresh=1e-7)$beta)
    ## scored X
    scoredX[,i] <- as.vector(centered.X %*% beta[,i])
    
    ## c. add new theta
    new.theta <- (diag(K) - as.matrix(Q[,i])%*%t(as.matrix(Q[,i]))%*%D) %*% solve(D)%*%t(Y)%*%X%*%beta[,i]
    new.norm <- as.numeric(t(as.matrix(new.theta))%*%D%*%as.matrix(new.theta))
    new.norm.theta <- new.theta * sqrt(1/new.norm)
    Q <- cbind(Q, new.norm.theta)
    ## are the set of scores mutually orthogonal?
    ## pracma::dot(Q[,2],Q[,3])
  }
  
  
  
  
  scoredX <- matrix(rep(0, nrow(Y)*(K-1)), ncol=K-1)
  
  ## scoredX with using trained beta 
  for(k in 1:(K-1)){
    scoredX[,k] <- as.vector(centered.X %*% beta[,k])        
  }
  
  ## LDA with two discriminant directions
  mu_k <- matrix(0, nrow = K-1, ncol = K) # centroid vectors
  pi_k <- vector(length = K) # priori probabilities
  sigma <- matrix(0, nrow = K-1, ncol = K-1) # covariance matrix
  
  for (l in 1:K){
    tmp <- matrix(scoredX[Y[,l]==1,],ncol=K-1)
    pi_k[l] <- nrow(tmp)/nrow(Y)
    mu_k[,l] <- colMeans(tmp)
    for(j in 1:nrow(tmp)){
      sigma <- sigma + (tmp[j,] - mu_k[,l]) %*% t(tmp[j,] - mu_k[,l])
    }
  }
  
  sigma <- sigma*(1/(nrow(Y)-K))
  inv.sigma <- solve(sigma)
  
  deltatrain <- matrix(0, ncol = K, nrow = nrow(Y))
  for (t in 1:K){
    deltatrain[,t] <- scoredX %*% inv.sigma %*% mu_k[,t]
    deltatrain[,t] <- deltatrain[,t] - 0.5*as.numeric(t(mu_k[,t]) %*% inv.sigma %*% mu_k[,t]) + log(pi_k[t])
  }
  
  
  return(model=list(beta=beta,inv.sigma =inv.sigma , mu_k= mu_k,pi_k=pi_k,mean_fac=scale.fac))
}
# opt.score(Y$train, X$train, R0, R1, lambda)