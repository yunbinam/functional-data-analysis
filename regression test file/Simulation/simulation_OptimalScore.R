# This simulation functions for optimal scoring methods
# ----------------------------------------------------------------------------

library(MASS)
simu_case1<-function(n=100,v1,v2,v3){
  # v1,v2,v3 eigen functions of laplace-baltrami opetrator
  # n is the number of samples in each class
  # train:valid:test=3:3:4
  
  Sigma<-diag(642)
    
  class1<-mvrnorm(n = n, v1, Sigma)
  class2<-mvrnorm(n = n, v2, Sigma)
  class3<-mvrnorm(n = n, v3, Sigma)

<<<<<<< HEAD

   X<-rbind(class1,class2,class3)
  label<-c(rep(0,n),rep(1,n),rep(2,n))
  indic1<-Reduce(rbind,lapply(1:n,function(x)matrix(c(1,0,0),ncol=3)))
  indic2<-Reduce(rbind,lapply(1:n,function(x)matrix(c(0,1,0),ncol=3)))
  indic3<-Reduce(rbind,lapply(1:n,function(x)matrix(c(0,0,1),ncol=3)))
  indic<-rbind(indic1,indic2,indic3)
  return(list(X=X,label=label,indic=indic))
}

data<-simu_case1(n=100,rep(1,642),rep(1,642),rep(1,642))


=======
  train<-rbind(class1[1:n1,],class2[1:n1,],class3[1:n1,])
  valid<-rbind(class1[(n1+1):n2,],class2[(n1+1):n2,],class3[(n1+1):n2,])
  test<-rbind(class1[(n2+1):n,],class2[(n2+1):n,],class3[(n2+1):n,])
  
  data<-list(train=train,valid=valid,test=test)
  
  return(data)
}

data<-simu_case1(n=100,v1,v2,v3)
>>>>>>> a63d91321cbb34fdce3c6ce00b870438c2ee59a6
