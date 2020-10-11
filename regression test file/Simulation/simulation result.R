# Arguments:
# sample:       generated samples, matrix A: [y X]
# R0:           mass matrix    
# R1:           stiffness matrix
# p:            proportion of training samples
# K:            the number of cross validation folds
# seq.l:        sequence of tuning parameter lambda
# random.folds: logical

oneSim <- function(sample, R0, R1,
                   p = 0.5, K = 5, seq.l = 10^(seq(-5,5,1))){
        y <- sample$y
        X <- sample$X
        
        n <- length(y)
        
        train <- sample(1:n, p * n)
        test <- -train
        
        X_train <- X[train,]
        X_test <- X[test,]
        y_train <- y[train]
        y_test <- y[test]
        
        cv.fit <- cv_spareg(y_train,X_train,R0,R1,K,seq.l)
        l.opt <- cv.fit$min_lambda
        
        coef <- spareg(y_train, X_train, R0, R1, l.opt)
        predict <- X_test%*%coef$beta+coef$intercept
        mse <- mean((y_test-predict)^2)
        
        return(mse)
}

        
        
