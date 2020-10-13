library(Matrix)
library(glmnet)
library(R.matlab)
library(ggplot2)
library(dplyr)
#path='\\\\fs2-vip/students/wenboz4/Desktop/high dimension project/DataAD/control group data/spareg'
#path='/Users/yunbi/Library/Mobile Documents/com~apple~CloudDocs/20_6_summer/Eardi/brain/regression test file'
#setwd(path)
source('spareg.R')
source('smooth_x.R')
source('smooth_recon_x.R')
source('CV_spareg.R')
source('cv_smooth_ols_x.R')
source('smooth_x_reg.R')
source('Simulation/sample generator.R')

# generate simulation samples
eigV <- readMat('Simulation/eigV.mat')$eigV
R0 <- readMM('R0_642.mtx')
R1 <- readMM('R1_642.mtx')

sample_noise <- mysamples(eigV, 10, 500, 4, 4, 1)
sample_smooth <- mysamples(eigV, 10, 500, 4, 0, 1)

#-----------------------------------------
#smoothing coefficients

mse1 <- replicate(100, oneSim(sample_noise, R0, R1, method="cv_spareg"))
mse2 <- replicate(100, oneSim(sample_noise, R0, R1, method="smooth_recon_x"))
mse3 <- replicate(100, oneSim(sample_noise, R0, R1, method="cv_smooth_ols_x"))

mse1_s <- replicate(100, oneSim(sample_smooth, R0, R1, method="cv_spareg"))
mse2_s <- replicate(100, oneSim(sample_smooth, R0, R1, method="smooth_recon_x"))
mse3_s <- replicate(100, oneSim(sample_smooth, R0, R1, method="cv_smooth_ols_x"))

rst <- data.frame("method"=c(rep("Smoothing beta", length(mse1)), 
                             rep("Smoothing x with recon", length(mse2)), 
                             rep("Smoothing x with cv", length(mse3))),
                  "sample"=c(rep("Noisy", length(mse1)+length(mse2)+length(mse3)),
                             rep("Smooth", length(mse1_s)+length(mse2_s)+length(mse3_s))),
                  "MSE"=c(mse1, mse2, mse3, mse1_s, mse2_s, mse3_s))

rst %>% group_by(sample) %>%
        ggplot(rst, aes(x=method, y=MSE, fill=method)) +
        geom_boxplot()


# compare beta* and beta
# integral t(beta*-beta)%*%R0%*%(beta*-beta) where beta: column vector
# box plot
# how much different box plots depending on the noise (on x and y) we put
# generalized eigenvalues 


