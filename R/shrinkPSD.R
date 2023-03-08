shrinkPSD <- function(R, 
                      n) {
    # Eigenvalues of the full covariance matrix
    lambdas = eigen(R)$values
    while ( min(lambdas) < 0 ) {
        R = 0.999*R	  # shrinkage
        diag(R) = 1			  # setting diagonal elements to 1
        lambdas = eigen(R)$values
    }
    C_XX_shr  =  R[1:n,  1:n]
    C_YY_shr  =  R[(n+1):dim(R)[1],  (n+1):dim(R)[1]]
    C_XY_shr  =  R[1:n,  (n+1):dim(R)[1]]
    return( list(C_XX_shr, C_YY_shr, C_XY_shr) )
} 