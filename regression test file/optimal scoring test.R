## beta from training data
beta <- cv_opt.score(Y$train, X$train, R0, R1, 5, 10^(seq(-5,5,1)))$beta

## test
X <- X$test
Y <- Y$test
K <- ncol(Y)

D <- 1/nrow(Y) * (t(Y)%*%Y) # K by K matrix

R0til.inv <- Matrix::Diagonal(x=1/Matrix::rowSums(R0))
L <- R1 %*% R0til.inv %*% R1
S <- Matrix::chol(L)
centered.X <- X - mean(X)

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
        tmp <- scoredX[Y[,l]==1,]
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

yhat <- apply(deltatrain, 1, which.max)
Y_label <- apply(Y, 1, which.max)
error <- mean(yhat!=Y_label)
