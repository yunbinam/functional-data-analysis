library(Matrix)
library(glmnet)
library(R.matlab)
library(ggplot2)
#path='\\\\fs2-vip/students/wenboz4/Desktop/high dimension project/DataAD/control group data/spareg'
#path='/Users/yunbi/Library/Mobile Documents/com~apple~CloudDocs/20_6_summer/Eardi/brain/regression test file'
#setwd(path)
source('spareg.R')
source('CV_spareg.R')
source('cv_smooth_x.R')
source('smooth_x_reg.R')
source('Simulation/sample generator.R')

# generate simulation samples
eigV <- readMat('Simulation/eigV.mat')$eigV
R0 <- readMM('R0_642.mtx')
R1 <- readMM('R1_642.mtx')

sample <- mysamples(eigV, 15, 500, 4, 1, 9)

#-----------------------------------------
#smoothing coefficients


mse1 <- replicate(100, simulation(sample, R0, R1, random.folds=TRUE))
mse2 <- replicate(100, simulation(sample, R0, R1, random.folds=FALSE))
mse3 <- replicate(100, oneSim(sample, R0, R1))
rst <- data.frame("f"=c(rep("f1", length(mse1)), rep("f2", length(mse2)), rep("f3", length(mse3))),
                  "MSE"=c(mse1, mse2, mse3))

ggplot(rst, aes(x=f, y=MSE)) +
        geom_boxplot()


# compare beta* and beta
# integral t(beta*-beta)%*%R0%*%(beta*-beta) where beta: column vector
# box plot
# how much different box plots depending on the noise (on x and y) we put
# generalized eigenvalues 


