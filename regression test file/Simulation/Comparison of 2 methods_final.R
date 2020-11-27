library(R.matlab)
library(Matrix)

source('../opt.score_new1.R')
source('../cv_opt.score_new1.R')
source('../test and results.R')

eigV <- readMat('eigV.mat')$eigV
R0 <- readMM('../R0_642.mtx')
R1 <- readMM('../R1_642.mtx')

## sample generator
mysamples <- function(eigV, K=3, n=100){
        
        r <- nrow(eigV) 
        c <- ncol(eigV)
        X <- matrix(0, nrow=(K*n), ncol=r) # n samples by r nodes 
        
        for(k in 1:K){
                for(i in 1:n){
                        # sample different bases v_1, v_2, v_3 for each individual 
                        base <- eigV[, 1:3]
                        # sample 3 by 1 mean 0 coefficients of bases
                        coef <- matrix(c(rnorm(1,0,5), rnorm(1,0,3), rnorm(1,0,1)), nrow=3)
                        # add different mean vector to each class
                        mean_k <- 0.3*eigV[,(3+k)]
                        # the distribution of coefficients are same in the same class 
                        x <- base%*%coef + rnorm(642, 0, 0.1) + mean_k
                        X[(k-1)*n+i,] <- x
                }
        }
        
        y <- Reduce(rbind, sapply(1:K, function(x) rep(x,n)))
        data <- cbind(y=y, X)
        new_data <- data[sample(nrow(data)),]
        
        y <- new_data[,1]
        X <- new_data[,2:ncol(new_data)]
        
        y_indice<-Reduce(rbind,lapply(y,function(x) {
                y_i <- rep(0,K)
                y_i[x]=1
                return(y_i)
        }))
        
        df <- data.frame(X=X, y_indice=y_indice, y=y)
        
        return(df)
}

## pre.smooth
pre.smooth <- function(ytrain, ytest, Xtrain, Xtest, R0, R1, lambdas){
        
        n_train <- nrow(Xtrain)
        n_test <- nrow(Xtest)
        p <- ncol(Xtrain)
        smooth.X_train <- matrix(0, nrow=n_train, ncol=p)
        smooth.X_test <- matrix(0, nrow=n_test, ncol=p)
        identity <- Diagonal(x=rep(1,p))
        
        R0til.inv <- Matrix::Diagonal(x=1/Matrix::rowSums(R0))
        S <- R1 %*% sqrt(R0til.inv)
        
        max <- 0
        select_lam <- 0
        lda_model <- NULL
        for(j in 1:length(lambdas)){
                lambda <- lambdas[j]
                Xstar <- rbind(identity, sqrt(lambda)*t(S))
                for(i in 1:n_train){
                        x_obs <- c(Xtrain[i,], rep(0, p))
                        smooth.X_train[i,] <- as.vector(glmnet(Xstar, x_obs, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta)
                }
                
                
                df.train <- data.frame(cbind(smooth.X_train, y=ytrain))
                
                fit <- lda(y~., df.train)
                fitted.y <- predict(fit, df.train)$class
                accuracy <- mean(ytrain == fitted.y)
                if(accuracy>max) {
                        max <- accuracy
                        select_lam <- lambda
                        lda_model <- fit
                }
        }
        
        Xstar <- rbind(identity, sqrt(select_lam)*t(S))
        for(i in 1:n_test){
                x_obs <- c(Xtest[i,], rep(0, p))
                smooth.X_test[i,] <- as.vector(glmnet(Xstar, x_obs, lambda = 0, alpha=0,intercept = FALSE, standardize = FALSE, thresh = 1e-7)$beta)
        }
        
        df.test <- data.frame(cbind(smooth.X_test, y=ytest))
        
        fitted.y <- predict(lda_model, df.test)$class
        accuracy <- mean(ytest == fitted.y)
        
        return(list(yhat = fitted.y, accuracy = accuracy))
}

my.lambda.seq <- 10^seq(-5,5,1)

result1.100 <- rep(0,100)

for(i in 1:100){
        
        samples <- mysamples(eigV,3,100)
        samples <- samples[order(samples$y),]
        
        X.train <- as.matrix(samples[,1:642])
        Y.train <- as.matrix(samples[,c("y_indice.1","y_indice.2","y_indice.3")])
        y.train <- as.numeric(samples[,"y"])
        
        samples.test <- mysamples(eigV,3,100)
        samples.test <- samples.test[order(samples.test$y),]
        
        X.test <- as.matrix(samples.test[,1:642])
        Y.test <- as.matrix(samples.test[,c("y_indice.1","y_indice.2","y_indice.3")])
        y.test <- as.numeric(samples.test[,"y"])
        
        fit.cv.opt <- cv_opt.score(Y.train, X.train, R0, R1, 5, my.lambda.seq)
        model.sm.flda <- opt.score(Y.train, X.train, R0, R1, fit.cv.opt$min_lambda)
        result.sm.flda <- test_acc(model.sm.flda, X.test, Y.test)
        result1.100[i] <- result.sm.flda$accuracy
}

boxplot(result1.100)

result2.100 <- rep(0, 100)

for(i in 1:100){
        
        samples <- mysamples(eigV,3,100)
        samples <- samples[order(samples$y),]
        
        X.train <- as.matrix(samples[,1:642])
        Y.train <- as.matrix(samples[,c("y_indice.1","y_indice.2","y_indice.3")])
        y.train <- as.numeric(samples[,"y"])
        
        samples.test <- mysamples(eigV,3,100)
        samples.test <- samples.test[order(samples.test$y),]
        
        X.test <- as.matrix(samples.test[,1:642])
        Y.test <- as.matrix(samples.test[,c("y_indice.1","y_indice.2","y_indice.3")])
        y.test <- as.numeric(samples.test[,"y"])
        
        result.presmooth <- pre.smooth(y.train, y.test, X.train, X.test, R0, R1, my.lambda.seq)
        result2.100[i] <- result.presmooth$accuracy

}

boxplot(result2.100)


result1.50 <- rep(0,50)

for(i in 1:50){
        
        samples <- mysamples(eigV,3,50)
        samples <- samples[order(samples$y),]
        
        X.train <- as.matrix(samples[,1:642])
        Y.train <- as.matrix(samples[,c("y_indice.1","y_indice.2","y_indice.3")])
        y.train <- as.numeric(samples[,"y"])
        
        samples.test <- mysamples(eigV,3,50)
        samples.test <- samples.test[order(samples.test$y),]
        
        X.test <- as.matrix(samples.test[,1:642])
        Y.test <- as.matrix(samples.test[,c("y_indice.1","y_indice.2","y_indice.3")])
        y.test <- as.numeric(samples.test[,"y"])
        
        fit.cv.opt <- cv_opt.score(Y.train, X.train, R0, R1, 5, my.lambda.seq)
        model.sm.flda <- opt.score(Y.train, X.train, R0, R1, fit.cv.opt$min_lambda)
        result.sm.flda <- test_acc(model.sm.flda, X.test, Y.test)
        result1.50[i] <- result.sm.flda$accuracy
}

boxplot(result1.50)

result2.50 <- rep(0, 50)

for(i in 1:50){
        
        samples <- mysamples(eigV,3,50)
        samples <- samples[order(samples$y),]
        
        X.train <- as.matrix(samples[,1:642])
        Y.train <- as.matrix(samples[,c("y_indice.1","y_indice.2","y_indice.3")])
        y.train <- as.numeric(samples[,"y"])
        
        samples.test <- mysamples(eigV,3,50)
        samples.test <- samples.test[order(samples.test$y),]
        
        X.test <- as.matrix(samples.test[,1:642])
        Y.test <- as.matrix(samples.test[,c("y_indice.1","y_indice.2","y_indice.3")])
        y.test <- as.numeric(samples.test[,"y"])
        
        result.presmooth <- pre.smooth(y.train, y.test, X.train, X.test, R0, R1, my.lambda.seq)
        result2.50[i] <- result.presmooth$accuracy
        
}

boxplot(result2.50)

df<- data.frame(accuracy=c(result1.50, result1.100, result2.50, result2.100), n=c(rep(150,50),rep(300,100),rep(150,50),rep(300,100)), method=c(rep("SM-FLDA",150), rep("Pre-smoothing",150)), method_n=c(rep("SM-FLDA, n=150", 50),rep("SM-FLDA, n=300",100),rep("Pre-smoothing, n=150", 50),rep("Pre-smoothing, n=300",100)))
