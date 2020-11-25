library(R.matlab)
library(Matrix)

source('sample generator.R')
source('../opt.score_new1.R')
source('../cv_opt.score_new1.R')
source('../test and results.R')
source('../presmoothing_lda.R')

eigV <- readMat('eigV.mat')$eigV
R0 <- readMM('../R0_642.mtx')
R1 <- readMM('../R1_642.mtx')

sample_c1 <- mysamples(eigV, R0, coef.mean=50, coef.var=1, noise_x=0.1)
sample_c2 <- mysamples(eigV, R0, coef.mean=30, coef.var=2, noise_x=0.1)
sample_c3 <- mysamples(eigV, R0, coef.mean=10, coef.var=3, noise_x=0.1)

X <- rbind(sample_c1, sample_c2, sample_c3)
Y <- rbind(matrix(rep(c(1,0,0),100), ncol=3, byrow=TRUE),
           matrix(rep(c(0,1,0),100), ncol=3, byrow=TRUE),
           matrix(rep(c(0,0,1),100), ncol=3, byrow=TRUE))
y <- c(rep(1,100), rep(2,100), rep(3,100))

X.train <- X[c(1:66,101:166,201:266),]
Y.train <- Y[c(1:66,101:166,201:266),]

X.test <- X[-c(1:66,101:166,201:266),]
Y.test <- Y[-c(1:66,101:166,201:266),]

my.lambda.seq <- 10^seq(-5,5,1)

fit.cv.opt <- cv_opt.score(Y.train, X.train, R0, R1, 5, my.lambda.seq)
model.sm.flda <- opt.score(Y.train, X.train, R0, R1, fit.cv.opt$min_lambda)
result.sm.flda <- test_acc(model.sm.flda, X.test, Y.test)

result.presmooth <- pre.smooth(y, X, R0, R1, 1)
