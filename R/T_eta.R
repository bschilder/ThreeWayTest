#' T_eta
#' 
#' Calculation of T_eta with eta fixed.
#' @param z_vector Column vectorized data matrix with rows represent 
#' phenotype and columns represent genotype.
#' @param cov_mat Estimated covariance matrix of z_vector.
#' @param eta Truncated value.
#' @returns A numeric value represent calculated statistic T_eta.
#' @importFrom MASS ginv
#' @export
#' @examples
#' z_vector<-MASS::mvrnorm(1,mu=rep(0,9),Sigma = diag(nrow = 9, ncol = 9))
#' ThreeWayTest::T_eta(z_vector=z_vector, 
#' cov_mat=diag(nrow = 9, ncol = 9), eta=0.5)
T_eta<-function(z_vector,
                cov_mat,
                eta){
    index<-ceiling(length(z_vector)*eta)
    z_vector_order<-order(abs(z_vector),decreasing = TRUE)
    ordered_z_vector<-z_vector[z_vector_order]
    ordered_cov_mat<-cov_mat[z_vector_order,z_vector_order]
    select_vector<-ordered_z_vector[1:index]
    select_cov_mat<-ordered_cov_mat[(1:index),(1:index)]
    stat<-t(select_vector)%*%MASS::ginv(select_cov_mat)%*%select_vector
    stat<-as.numeric(stat)
    return(stat)
}