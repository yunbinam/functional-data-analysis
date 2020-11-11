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
eigV <- R.matlab::readMat('Simulation/eigV.mat')$eigV
v <- eigV[,sample(1:100,3)]
X <- simu_case1(n=100, v1=v[,1], v2=v[,2], v3=v[,3])
Y <- list(train=rbind(matrix(rep(c(1,0,0),30), ncol=3, byrow=TRUE),
                      matrix(rep(c(0,1,0),30), ncol=3, byrow=TRUE),
                      matrix(rep(c(0,0,1),30), ncol=3, byrow=TRUE)),
          test=rbind(matrix(rep(c(1,0,0),40), ncol=3, byrow=TRUE),
                      matrix(rep(c(0,1,0),40), ncol=3, byrow=TRUE),
                      matrix(rep(c(0,0,1),40), ncol=3, byrow=TRUE)))

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
                centered.mapped.Y <- mapped.Y - mean(mapped.Y)
                Y.star <- c(centered.mapped.Y, rep(0, s))
                
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
        
        return(beta)
}

# opt.score(Y$train, X$train, R0, R1, lambda)