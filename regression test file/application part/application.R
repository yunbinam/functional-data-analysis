library(Rvcg)
library(rgl)
library(shapes)
library(Matrix)
library(dplyr)
library(igraph)
library(PEIP)
library(glmnet)
library(R.matlab)

#------------------------------------------------------------------------------
# axuliary functions
test_acc<-function(model,x_test,y_test){
  K=ncol(y_test)
  beta=model$beta
  mu_k <- model$mu_k # centroid vectors
  pi_k <- model$pi_k # priori probabilities
  inv.sigma <- model$inv.sigma # covariance matrix
  
  center.test<-sapply(1:length(model$mean_fac),function(x) x_test[,x]-model$mean_fac[x])
  scoredX<-center.test%*%beta
  
  deltatrain <- matrix(0, ncol = K, nrow = nrow(y_test))
  for (t in 1:K){
    deltatrain[,t] <- scoredX %*% inv.sigma %*% mu_k[,t]
    deltatrain[,t] <- deltatrain[,t] - 0.5*as.numeric(t(mu_k[,t]) %*% inv.sigma %*% mu_k[,t]) + log(pi_k[t])
  }
  
  yhat <- apply(deltatrain, 1, which.max)
  Y_label <- apply(y_test, 1, which.max)
  accuracy<- mean(yhat==Y_label)
  list(accuracy=accuracy,yhat=yhat)
}

get_indice<-function(y){
  Reduce(rbind,lapply(1:length(y),function(i){
    if(y[i]==0) return(matrix(c(1,0),ncol=2))
    else return(matrix(c(0,1),ncol=2))
  }))
}
#---------------------------------------------------------------------------------------

# set path
setwd("/Users/wenbozhang/Desktop/high dimensional/application part")


p=32492
shape_files<- grep('L\\.midthickness\\.32k_fs_LR\\.ply$', list.files("files1"), value=TRUE)
thickness_files<-grep('L\\.thickness\\.32k_fs_LR\\.func\\.1D', list.files("files1"), value=TRUE)
id<- gsub(".L.+", "", shape_files)
for(i in 1:length(id)){
  if(substr(id[i],0,1)=="P"){
    number<-substr(id[i],3,6)
    id[i]<-paste0("P",number)
  }
}

## clean the grouo information
group_info <- read.csv('All_grp_ATN.csv', header=TRUE)
group_info<-group_info[group_info$Grp!="NAD",]
group_info$Grp<-as.character(group_info$Grp)

select_sample<-id[id%in%group_info$Subject]
label<-rep(0,length(select_sample))

## two reptition samples 165-2=-163
for(i in 1:length(label)){
  class<-group_info[group_info$Subject==select_sample[i],][1,]$Grp
  label[i]<-ifelse(class=="C",0,1)
}

path<-thickness_files[id%in%select_sample]
fdata = t(sapply(path, function(x) as.matrix(read.csv(paste0('files1/',x), header = F))))

## shuffle data
data<-cbind(label,fdata)

#divide train and test dataset
data_samples<-function(data){
  data2<-data[sample(nrow(data)),]
  x_data<-data2[,2:32493]
  y<-data2[,1]
  
  y_indice<-get_indice(y)
  
  train_size=114
  x_train<-x_data[1:train_size,]
  x_test<-x_data[(train_size+1):nrow(x_data),]
  y_train<-y_indice[1:train_size,]
  y_test<-y_indice[(train_size+1):nrow(y_indice),]
  
  return(list(x_train=x_train,x_test=x_test,y_train=y_train,y_test=y_test))
}

# read r0 r1
R0=readMM('R0.mtx')
R1=readMM('R1.mtx')

#R0=readMM('R0_new.mtx')
#R1=readMM('R1_new.mtx')

#R0_new_scale=readMM('R0_new_scale.mtx')
#R1_new_scale=readMM('R1_new_scale.mtx')



samples<-data_samples(data)
x_train<-samples$x_train
x_test<-samples$x_test
y_train<-samples$y_train
y_test<-samples$y_test

cv_model<-cv_opt.score(y_train, x_train, R0, R1, 5,10^seq(-5,5,0.5))
model1<-opt.score(y_train, x_train, R0, R1, 0.1)
result1<-test_acc(model1,x_test,y_test)
result1$accuracy
result1$yhat


# cvfit1<-cv.glmnet(x_train, apply(y_train, 1, which.max), nfolds=5,alpha=0,lambda=10^seq(-5,5,0.5),family = "binomial", type.measure = "class")
# pre<-predict(cvfit, newx = x_test, s = "lambda.min", type = "class")
# mean(as.numeric(pre==apply(y_test, 1, which.max)))
# 
# cvfit2<-cv.glmnet(x_train, apply(y_train, 1, which.max), nfolds=5,alpha=0.5,lambda=10^seq(-5,5,0.5),family = "binomial", type.measure = "class")
# pre2<-predict(cvfit2, newx = x_test, s = "lambda.min", type = "class")
# mean(as.numeric(pre2==apply(y_test, 1, which.max)))
# 
# cvfit3<-cv.glmnet(x_train, apply(y_train, 1, which.max), nfolds=5,alpha=1,lambda=10^seq(-5,5,0.5),family = "binomial", type.measure = "class")
# pre3<-predict(cvfit3, newx = x_test, s = "lambda.min", type = "class")
# mean(as.numeric(pre3==apply(y_test, 1, which.max)))




