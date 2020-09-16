library(Matrix)


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

# modeling
start.time <- Sys.time()
coef=spareg(y_train,X_train,R0,R1,0.01)
end.time <- Sys.time()
print(end.time-start.time)

#--------------------------------------
# performance
predict=X_test%*%coef$beta+coef$intercept
mean((y_test-predict)^2)


## -------------------------------
# compare with logistic regression with glmnet
model1=glmnet(X_train, y_train, lambda = 0, alpha=0,intercept = TRUE, standardize = TRUE, thresh = 1e-7)
predict=X_test%*%model1$beta+model1$a0
mean((y_test-predict)^2)




##-----------
# cross-validation
reg=cv_spareg(Y, X, R0,R1, 2, c(0.01,0.02))






