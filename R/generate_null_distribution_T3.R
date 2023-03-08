#' generate_null_distribution_T3
#' 
#' This function generate the null distribution of T_eta,
#'  should be run before calculating TWT.
#' Please see Setp 2 in our manuscript for more information.
#' @param m Dimension.
#' @param n Number of the generated sample size.  
#' @param cov_mat Covariance matrix for generating random samples.
#' @param cutoff_value Truncated value, should be in 0 to 1 and must contain 1.
#' @returns Null distribution matrix with number of rows equal n and number of 
#' columns equal length of cutoff value minus 1.
#' @importFrom MASS mvrnorm
#' @export
#' @examples 
#' ThreeWayTest::generate_null_distribution_T3(m=6,n=1000,
#' cov_mat=diag(nrow = 6, ncol = 6), cutoff_value=c(0.2,0.4,0.6,0.8,1))
generate_null_distribution_T3<-function(m,
                                        n,
                                        cov_mat,
                                        cutoff_value){
    true_mean<-rep(0,m)
    bootstrap_data<-MASS::mvrnorm(n,true_mean,cov_mat)
    null_distribution_matrix<-matrix(0,nrow = n, ncol = length(cutoff_value)-1)
    for (i in 1:(length(cutoff_value)-1)) {
        for (j in 1:n) {
            null_distribution_matrix[j,i]<-T_eta(z_vector = bootstrap_data[j,],
                                                 cov_mat,cutoff_value[i])
        }
    }
    return(null_distribution_matrix)
}