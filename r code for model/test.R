library(Rvcg)
library(rgl)
library(shapes)
library(Matrix)
library(dplyr)
library(glmnet)
library(R.matlab)
library(glmnet)
library(fdaPDE)
library(MASS)



#------------------------------------------------------------------------------------
# simulation function
mysamples <- function(eigV,K=3, n = 100){
  r <- nrow(eigV)
  c <- ncol(eigV)
  X <- matrix(0, nrow=(K*n), ncol=r) # n samples by r nodes 
  for(k in 1:K){
    for(i in 1:n){
      # sample different bases v_1, v_2, v_3 for each individual 
      base <- eigV[, 1:3]
      coef <- matrix(c(rnorm(1, 0,5),rnorm(1, 0,3),rnorm(1, 0,1)), nrow=3)
      mean_k <- 0.3*eigV[,(3+k)]
      # the distribution of coefficients are same in the same class 
      x <- base%*%coef + rnorm(642, 0, 0.1)+mean_k
      X[(k-1)*n+i,] <- x
    }
  }
  y<-Reduce(rbind,sapply(1:K,function(x)rep(x,n)))
  data <- cbind(y=y,X)
  new_data<-data[sample(nrow(data)),]
  
  y<-new_data[,1]
  X<-new_data[,2:ncol(new_data)]
  y_indice<-Reduce(rbind,lapply(y,function(x) {
    y_i <- rep(0,K)
    y_i[x]=1
    return(y_i)
  }))
  return(list(X=X,y_indice=y_indice))
}
#----------------------------------------------------------------------------

setwd('/Users/wenbozhang/Desktop/high dimensional/application part/r code for model')

# read r0 r1, eigv
R0_642=readMM('R0_642.mtx')
R1_642=readMM('R1_642.mtx')
eigv<-readMat("eigv.mat")$eigV

# FEM_basis
nodes_642<-readMat('manifold_642.mat')$manifold[1][[1]]
trian_642<-readMat('manifold_642.mat')$manifold[2][[1]]
mesh_642 = create.mesh.2.5D(nodes = nodes_642, triangles = trian_642)
FEMbasis_642 = create.FEM.basis(mesh_642)

set.seed(1)
# simulate samples
samples<-mysamples(eigv,3,100)
X<-samples$X
Y_indice<-samples$y_indice

x_train<-X[1:210,]
y_train<-Y_indice[1:210,]
x_test<-X[211:300,]
y_test<-Y_indice[211:300,]

cv_model<-cv.OptScore(x_train,y_train,R0_642, R1_642, 5,10^seq(-3,3,1))
model1<-OptScore(x_train, y_train, R0_642, R1_642, cv_model$min_lambda)
predict(model1,x_test,y_test)

# here the dimension d of FPCA is 10
cv_fpca<-cv.FPCA_LDA(x_train,apply(y_train, 1, which.max),FEMbasis_642,lambdas=10^seq(-3,3,1),kfolds=5,d=10)
model2<-FPCA_LDA(x_train,apply(y_train, 1, which.max),FEMbasis_642,cv_fpca$min_lambda,15)
predict(model2,x_test,apply(y_test, 1, which.max))

cv_pre<-cv.PreSmooth(x_train,apply(y_train, 1, which.max),R0_642, R1_642, 5,10^seq(-3,3,1))
model3<-PreSmooth(x_train,apply(y_train, 1, which.max),R0_642, R1_642,cv_pre$min_lam)
predict(model3,x_test,apply(y_test, 1, which.max),cv_pre$min_lam)


