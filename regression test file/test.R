library(Matrix)
library(glmnet)
library(R.matlab)
#path='\\\\fs2-vip/students/wenboz4/Desktop/high dimension project/DataAD/control group dataX/spareg'
setwd(path)
source('spareg.R')
source('cv_spareg.R')
source('cv_smooth_ols_x.R')
source('smooth_x.R')
source('smooth_recon_x.R')

# reading the data 
#data=readRDS('A.rds')
#R0=readMM('R0.mtx')
#R1=readMM('R1.mtx')

# reading the data 
data=readMat('A_samples_642.mat')
R0=readMM('R0_642.mtx')
R1=readMM('R1_642.mtx')


y=data$A[,1]
X=data$A[,2:ncol(data$A)]

n=length(y)

set.seed(10)
train=sample(1:length(y), 0.8 * length(y))
test=-train

X_train=X[train,]
X_test=X[test,]
y_train=y[train]
y_test=y[test]



## 5 functions used here: 1,spareg 2.cv_spareg 3.smooth_x 4.cv_smooth_ols_x 5.smooth_recon_x:
####### 1 and 2 do a regression with smoothing the coefficients beta; 
####### 3 to 5 smooth covariates x first and then do a regression
# Notice:
# 4 and 5 can be seen as two different ways to select lambda, the difference is 4 minimizes cv mse and 5 
# minimizes reconstruction error of X

#-----------------------------------------
#1.smoothing coefficients with fixed lambda
coef1=spareg(y_train,X_train,R0,R1,0.01)
predict1=X_test%*%coef1$beta+coef1$intercept
mse1=mean((y_test-predict1)^2)

##-----------------------------
#2.smoothing coefficients with cross-validation
reg=cv_spareg(y_train,X_train,R0,R1,5,10^(seq(-3,3,1)))

# using selected lambda to train the model
coef2=spareg(y_train,X_train,R0,R1,reg$min_lambda)

predict2=X_test%*%coef2$beta+coef2$intercept
mse2=mean((y_test-predict2)^2)



##-----------------------------
#3.smoothing covariates with fixed lambda
coef3=smooth_x(y_train,X_train,R0,R1,1)

# using selected lambda to smooth test dataset first
smooth_X_test=smooth_x(y_test,X_test,R0,R1,1,'test')

predict3=smooth_X_test%*%coef3$beta+coef3$intercept
mse3=mean((y_test-predict3)^2)


##-----------------------------
#4.smoothing covariates with cross validation which minimize least squares
model1=cv_smooth_ols_x(y_train,X_train, R0,R1, 5, 10^(seq(-3,3,1)))
coef4=model1$coef

# using selected lambda to smooth test dataset first
smooth_X_test=smooth_x(y_test,X_test,R0,R1,model$min_lambda,'test')
predict4=smooth_X_test%*%coef4$beta+coef4$intercept
mse4=mean((y_test-predict4)^2)


##-----------------------------
#5.smoothing covariates with the lambda minimizing reconstruction error
model2=smooth_recon_x(y_train,X_train, R0,R1, 10^(seq(-3,3,1)))
coef5=model2$coef

# using selected lambda to smooth test dataset first
smooth_X_test=smooth_x(y_test,X_test,R0,R1,model2$min_lambda,'test')
predict5=smooth_X_test%*%coef5$beta+coef5$intercept
mse5=mean((y_test-predict5)^2)
