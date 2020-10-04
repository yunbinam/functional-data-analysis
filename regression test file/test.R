library(Matrix)
library(glmnet)
#path='\\\\fs2-vip/students/wenboz4/Desktop/high dimension project/DataAD/control group data/spareg'
#setwd(path)
source('spareg.R')
source('CV_spareg.R')
source('cv_smooth_x.R')
source('smooth_x_reg.R')

# reading the data 
data=readRDS('A.rds')
R0=readMM('R0.mtx')
R1=readMM('R1.mtx')

y=data[,1]
X=data[,2:32493]

X_train=X[1:70,]
X_test=X[70:87,]
y_train=y[1:70]
y_test=y[70:87]

#-----------------------------------------
#smoothing coefficients
coef=spareg(y_train,X_train,R0,R1,0.01)
predict=X_test%*%coef$beta+coef$intercept
mean((y_test-predict)^2)


# cross-validation
reg=cv_spareg(y_train, X_train, R0,R1, 5,10^(seq(-3,3,1)))
coef=spareg(y_train,X_train,R0,R1,reg$min_lambda)
predict=X_test%*%coef$beta+coef$intercept
mean((y_test-predict)^2)



##-----------------------------
# smoothing covariates
coef=smooth_x_reg(y[1:5],X[1:5,],R0,R1,1)

#cv
model=cv_smooth_x(y, X, R0,R1, 5, c(0,0.1))