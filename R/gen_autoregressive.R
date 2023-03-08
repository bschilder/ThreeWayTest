#' gen_autoregressive
#'
#' This is the function for generating an autoregressive covariance matrix.
#' @param m Dimension of the matrix.
#' @param rho Correlation coefficient.
#' @returns An autoregressive covariance matrix.
#' 
#' @export
#' @examples
#' cov_mat <- gen_autoregressive(m = 6, rho = 0.5)
gen_autoregressive<-function(m,
                             rho){
    cov_mat<-matrix(0,m,m)
    for (i in 1:m){
        for (j in 1:m){
            cov_mat[i,j]=rho^(abs(i-j))
        }
    }
    return(cov_mat)
}