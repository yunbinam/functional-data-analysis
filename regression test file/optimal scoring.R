# Arguments:
m# X:          n by s design matrix (covariates)
# R0:         s by s symmetric mass matrix
# R1:         s by s symmetric stiffness matrix
# lambda:     tuning parameter of smoothness penalty
# ------------------------------------------------------------------------
# Outputs:
# 

## simple example
Y <- matrix(c(0,1,0,
              1,0,0,
              1,0,0,
              0,0,1,
              0,0,1,
              1,0,0,
              0,1,0), ncol=3, byrow=TRUE)

X <- matrix(runif(7*642), ncol=642)

R0 <- Matrix::readMM('R0_642.mtx')
R1 <- Matrix::readMM('R1_642.mtx')
lambda <- 1

opt.score <- function(Y, X, R0, R1, lambda){
        
        Y <- as.matrix(Y)
        X <- as.matrix(X)
        
        n <- nrow(Y)
        K <- ncol(Y)
        s <- ncol(X)
        
        ## 1. prepare matrices
        D <- 1/n * (t(Y)%*%Y) # K by K matrix
        
        R0til.inv <- Matrix::Diagonal(x=1/Matrix::rowSums(R0))
        L <- R1 %*% R0til.inv %*% R1
        S <- Matrix::chol(L)
        centered.X <- X - mean(X)
        X.star <- rbind(centered.X, sqrt(lambda)*t(as.matrix(S)))
        
        ## 3. pre-specify theta1
        Q <- as.matrix(rep(1, ncol(Y)))
        predicted.scores <- matrix(rep(0, n*(K-1)), ncol=K-1)
        
        ## 4. loop
        for(i in 1:(K-1)){
                ## a. compute theta
                random.theta <- runif(K)
                theta <- (diag(K) - as.matrix(Q[,i])%*%t(as.matrix(Q[,i]))%*%D) %*% random.theta
                norm <- as.numeric(t(theta)%*%D%*%theta)
                norm.theta <- theta * sqrt(1/norm)
                
                ## b. solve beta
                mapped.Y <- Y %*% norm.theta
                centered.mapped.Y <- mapped.Y - mean(mapped.Y)
                Y.star <- c(centered.mapped.Y, rep(0, s))
                
                beta <- glmnet::glmnet(X.star, Y.star, 
                                       lambda=0, alpha=0, intercept=FALSE, 
                                       standardize=FALSE, thresh=1e-7)$beta
                ## predicted scores
                predicted.scores[,i] <- as.vector(centered.X %*% beta)
                
                ## c. add new theta
                new.theta <- (diag(K) - as.matrix(Q[,i])%*%t(as.matrix(Q[,i]))%*%D) %*% solve(D)%*%t(Y)%*%X%*%beta
                new.norm <- as.numeric(t(as.matrix(new.theta))%*%D%*%as.matrix(new.theta))
                new.norm.theta <- new.theta * sqrt(1/new.norm)
                Q <- cbind(Q, new.norm.theta)
                ## are the set of scores mutually orthogonal?
                ## pracma::dot(Q[,2],Q[,3])
        }
        
        ## LDA with the first (one) discriminant direction
        mu_k <- matrix(0, nrow = s, ncol = K) # centroid vectors
        pi_k <- vector(length = K)
        sigma <- matrix(0, nrow = s, ncol = s) # covariance matrix
        
        for (i in 1:K){
                tmp <- X[Y[,i]==1,]
                pi_k[i] <- nrow(tmp)/n
                mu_k[,i] <- colMeans(tmp)
                for(j in 1:nrow(tmp)){
                        sigma <- sigma + (tmp[j,] - mu_k[,i]) %*% t(tmp[j,] - mu_k[,i])
                }
        }
        
        sigma <- sigma*(1/(n-K))
        inv.sigma <- solve(sigma) # error.. not invertible
        
        ## How do we use score functions?
        
}

