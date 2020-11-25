# Arguments:
# eigV:       r by c matrix, c eigenfunctions of the Laplace-Beltrami operator
# m:          the number of base eigenfunctions using in linear combination to generate x_ii
# n:          the number of generating samples
# coef.mean:  mean of coefficients in linear combination
# coef.var:   variance of coefficients in linear combination
# noise_x:    varaince of random noise on x
# ------------------------------------------------------------------------
# Outputs:
# A:          n by (r+1) matrix [y X]
# y:          n by 1 matrix 
# X:          n by r matrix 

mysamples <- function(eigV, R0, m = 3, n = 100, coef.mean, coef.var, noise_x){
        
        r <- nrow(eigV)
        c <- ncol(eigV)
        
        X <- matrix(0, nrow=n, ncol=r) # n samples by r nodes
        
        for(i in 1:n){
                
                base_index <- sample(1:c,m) 
                # sample different bases v_1, v_2, v_3 for each individual
                base <- eigV[, base_index]
                
                coef <- matrix(rnorm(m, coef.mean, coef.var), nrow=m) 
                # the distribution of coefficients are same in the same class
                x <- base%*%coef + rnorm(r, 0, noise_x)
                X[i,] <- x
                
        }
        
        return(X)
}