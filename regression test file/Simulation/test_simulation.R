library(Matrix)
library(glmnet)
library(R.matlab)
#path='\\\\fs2-vip/students/wenboz4/Desktop/high dimension project/DataAD/control group data/spareg'
#path='/Users/yunbi/Library/Mobile Documents/com~apple~CloudDocs/20_6_summer/Eardi/brain/regression test file'
#setwd(path)
source('spareg.R')
source('CV_spareg.R')
source('cv_smooth_x.R')
source('smooth_x_reg.R')

# reading the data 
data=readMat('Simulation/A_samples.mat')
R0=readMM('R0.mtx')
R1=readMM('R1.mtx')

y=data$A[,1]
X=data$A[,2:32493]

n=length(y)
train=sample(1:length(y), 0.8 * length(y))
test=-train

X_train=X[train,]
X_test=X[test,]
y_train=y[train]
y_test=y[test]

#-----------------------------------------
#smoothing coefficients
coef=spareg(y_train,X_train,R0,R1,0.01)
predict=X_test%*%coef$beta+coef$intercept
mse1=mean((y_test-predict)^2)


# cross-validation
reg=cv_spareg(y_train, X_train, R0,R1, 5,10^(seq(-3,3,1)))
coef=spareg(y_train,X_train,R0,R1,reg$min_lambda)
predict=X_test%*%coef$beta+coef$intercept
mse2=mean((y_test-predict)^2)



##-----------------------------
# smoothing covariates
# coef=smooth_x_reg(y[1:5],X[1:5,],R0,R1,1)

#cv
# model=cv_smooth_x(y, X, R0,R1, 5, c(0,0.1))