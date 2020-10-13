# Arguments:
# sample:       generated samples, matrix A: [y X]
# R0:           mass matrix    
# R1:           stiffness matrix
# method:       one of comparing methods: cv_spareg, smooth_recon_x, cv_smooth_ols_x
# p:            proportion of training samples
# K:            the number of cross validation folds
# seq.l:        sequence of tuning parameter lambda
# random.folds: logical

oneSim <- function(sample, R0, R1, method = "cv_spareg",
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
        
        if(method=="cv_spareg"){
                spareg.fit <- cv_spareg(y_train,X_train,R0,R1,K,seq.l)
                coef <- spareg.fit$coef
        }
        
        if(method=="smooth_recon_x"){
                recon.fit <- smooth_recon_x(y_train,X_train,R0,R1,seq.l)
                coef <- recon.fit$coef
        }
        
        if(method=="cv_smooth_ols_x"){
                ols.fit <- cv_smooth_ols_x(y_train,X_train,R0,R1,K,seq.l)
                coef <- ols.fit$coef
        }
        
        predict <- X_test%*%coef$beta+coef$intercept
        mse <- mean((y_test-predict)^2)
        
        return(mse)
}

        
        
