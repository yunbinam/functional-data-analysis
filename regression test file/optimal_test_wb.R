eigV <- R.matlab::readMat('eigV.mat')$eigV
v <- eigV[,sample(1:100,3)]

data<- simu_case1(n=100, v1=100*v[,1], v2=100*v[,2], v3=100*v[,3])
R0 <- Matrix::readMM('R0_642.mtx')
R1 <- Matrix::readMM('R1_642.mtx')

X<-data$X
indic<-data$indic
X_train<-X[c(1:30,101:130,201:230),]
Y_train<-indic[c(1:30,101:130,201:230),]
X_test<-X[-c(1:30,101:130,201:230),]
Y_test<-indic[-c(1:30,101:130,201:230),]

model<-opt.score(Y_train, X_train, R0, R1, 1)

K=3
beta=model$beta
mu_k <- model$mu_k # centroid vectors
pi_k <- model$pi_k # priori probabilities
inv.sigma <- model$inv.sigma # covariance matrix

scoredX<-cbind(X_test%*%beta[,1],X_test%*%beta[,2])



deltatrain <- matrix(0, ncol = K, nrow = nrow(Y_test))
for (t in 1:K){
  deltatrain[,t] <- scoredX %*% inv.sigma %*% mu_k[,t]
  deltatrain[,t] <- deltatrain[,t] - 0.5*as.numeric(t(mu_k[,t]) %*% inv.sigma %*% mu_k[,t]) + log(pi_k[t])
}

yhat <- apply(deltatrain, 1, which.max)
Y_label <- apply(Y_test, 1, which.max)
accuracy<- mean(yhat==Y_label)

plot_data<-data.frame(scoredX,label=Y_label)

library(ggplot2)
ggplot(plot_data, aes(x=X1, y=X2,color=label)) + 
  geom_point(size=3)
