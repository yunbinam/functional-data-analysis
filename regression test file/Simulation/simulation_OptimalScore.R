# This simulation functions for optimal scoring methods
# ----------------------------------------------------------------------------

library(MASS)
simu_case1<-function(n=100,v1,v2,v3){
  # v1,v2,v3 eigen functions of laplace-baltrami opetrator
  # n is the number of samples in each class

  
  Sigma<-diag(642)
  n1<-3*(n/10)
  n2<-2*n1
    
  class1<-mvrnorm(n = n, v1, Sigma)
  class2<-mvrnorm(n = n, v2, Sigma)
  class3<-mvrnorm(n = n, v3, Sigma)


  X<-rbind(class1,class2,class3)
  label<-c(rep(0,30),rep(1,30),rep(2,30))
  return(list(X=X,label=label))
}

data<-simu_case1(n=100,v1,v2,v3)


