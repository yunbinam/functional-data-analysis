# Arguments:
# eigV:       r by c matrix, c eigenfunctions of the Laplace-Beltrami operator
# m:          the number of base eigenfunctions using in linear combination to generate x_ii
# n:          the number of generating samples
# v:          variance of coefficients in linear combination
# noise_x:    varaince of random noise on x
# noise_y:    varaince of random noise on y
# ------------------------------------------------------------------------
# Outputs:
# A:          n by (r+1) matrix [y X]
# y:          n by 1 matrix 
# X:          n by r matrix 

mysamples <- function(eigV, m, n, v, noise_x, noise_y){
        
        r <- nrow(eigV)
        c <- ncol(eigV)
        
        # pick an eigenfunction as beta
        beta <- matrix(eigV[, sample(1:c,1)], ncol=1) # r by 1 vector
        X <- matrix(0, nrow=n, ncol=r) # n samples by r nodes
        
        for(i in 1:n){
                
                base_index <- sample(1:c,m)
                base <- eigV[, base_index]
                coef <- matrix(rnorm(m, 0, v), nrow=m)
                x <- base%*%coef + rnorm(r, 0, noise_x)
                X[i,] <- x
                
        }
        
        y <- X %*% beta + rnorm(n, 0, noise_y) # n by 1 vector
        
        A <- cbind(y, X)
        
        return(list(A=A, y=y, X=X))
}

