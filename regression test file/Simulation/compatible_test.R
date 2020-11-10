library(Matrix)
library(R.matlab)
library(ggplot2)

source('sample generator.R')

eigV <- readMat('eigV.mat')$eigV
R0 <- readMM('../R0_642.mtx')
R1 <- readMM('../R1_642.mtx')

## Generate our test sample
sample <- mysamples(eigV, 10, 500, 4, 1, 1)

## Calculate the area of manifold to normalize the sum
## By normalizing with the area, the sum reflects how big the manifold is
one_vector <- rep(1, nrow(R0))
manifold_area <- t(one_vector) %*% R0 %*% one_vector

## Use normalized sum in case we do not fully observe x_i
normalized_sum <- (manifold_area/nrow(R0))[1,1] * (sample$X %*% sample$beta)
## Use integral in case we have fully observed functional x_i
integral <- as.vector(sample$X %*% R0 %*% sample$beta)

results <- data.frame(formulation=c(rep("sum",500),rep("integral",500)), 
                      results=c(normalized_sum,integral))

## Check if the two formulations are asymptotically equivalent
## even though it does not compare whether the results are compatible at each point
## we can use it for the purpose of checking the range/scale of results from each formulation
ggplot(results, aes(x=formulation, y=results, color=results))+
        geom_boxplot()
