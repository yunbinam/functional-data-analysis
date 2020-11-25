# Arguments:
# y:          n by 1 class vector (response)
# X:          n by s design matrix (covariates)
# R0:         s by s symmetric mass matrix
# R1:         s by s symmetric stiffness matrix
# lambda:     tuning parameter of smoothness penalty
# ------------------------------------------------------------------------
# Outputs:
# smooth.X:   n by s smoothed X matrix

## simple example
library(glmnet)
library(MASS)

pre.smooth <- function(y, X, R0, R1, lambda){
        
        n <- nrow(X)
        p <- ncol(X)
        smooth.X <- matrix(0,nrow=n,ncol=p)
        identity <- Diagonal(x=rep(1,p))
        
        R0til.inv <- Matrix::Diagonal(x=1/Matrix::rowSums(R0))
        S <- R1 %*% sqrt(R0til.inv)
        
        Xstar <- rbind(identity, sqrt(lambda)*t(S))
        for(i in 1:n){
                x_obs <- c(X[i,], rep(0, p))
                smooth.X[i,] <- as.vector(glmnet(Xstar, x_obs, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta)
        }
        
        smooth.X.train <- smooth.X[c(1:66,101:166,201:266),]
        y.train <- y[c(1:66,101:166,201:266)]
        
        smooth.X.test <- smooth.X[-c(1:66,101:166,201:266),]
        y.test <- y[-c(1:66,101:166,201:266)]
        
        df.train <- data.frame(cbind(smooth.X.train, y=y.train))
        df.test <- data.frame(cbind(smooth.X.test, y=y.test))
        
        fit <- lda(y~., df.train, prior=c(1,1,1)/3)
        fitted.y <- predict(fit, df.test)$class
        accuracy <- mean(y.test == fitted.y)
        
        return(list(yhat = fitted.y, accuracy = accuracy))
}
        