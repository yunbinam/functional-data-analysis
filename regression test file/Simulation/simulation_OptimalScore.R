# This simulation functions for optimal scoring methods
# ----------------------------------------------------------------------------

library(MASS)
simu_case1<-function(n=100,v1,v2,v3){
  # v1,v2,v3 eigen functions of laplace-baltrami opetrator
  # n is the number of samples in each class
  # train:valid:test=3:3:4
  
  Sigma<-diag(642)
  n1<-3*(n/10)
  n2<-2*n1
    
  class1<-mvrnorm(n = n, v1, Sigma)
  class2<-mvrnorm(n = n, v2, Sigma)
  class3<-mvrnorm(n = n, v3, Sigma)

  train<-rbind(class1[1:n1,],class2[1:n1,],class3[1:n1,])
  valid<-rbind(class1[(n1+1):n2,],class2[(n1+1):n2,],class3[(n1+1):n2,])
  test<-rbind(class1[(n2+1):n,],class2[(n2+1):n,],class3[(n2+1):n,])
  
  data<-list(train=train,valid=valid,test=test)
  
  return(data)
}

data<-simu_case1(n=100,v1,v2,v3)
