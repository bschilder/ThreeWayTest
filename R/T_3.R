#' T_3
#' 
#' This is the function calculating T3.
#' @param z_vector Column vectorized data matrix with rows 
#' represent phenotype and columns represent genotype.
#' @param cov_mat Estimated covariance matrix of z_vector.
#' @param cutoff_value Vector of truncated value eta.
#' @param coefficient_matrix Matrix of alpha, beta and d.
#' @returns A numeric value represent p value of T3.
#' @importFrom MASS ginv
#' @export
#' @examples
#' z_vector<-MASS::mvrnorm(1,mu=rep(0,10),
#'                         Sigma = diag(nrow = 10, ncol = 10))
#' null_distribution<-generate_null_distribution_T3(
#'     m=10,
#'     n=1000,
#'     cov_mat=diag(nrow = 10,ncol= 10), 
#'     cutoff_value=c(0.2,0.4,0.6,0.8,1))
#' coefficient_matrix<- approximate_distribution_coefficient_estimate_T3(
#'         null_distribution_matrix = null_distribution)
#' pvalue <- T_3(z_vector=z_vector, 
#'               cov_mat=diag(nrow = 10, ncol = 10),
#'               cutoff_value=c(0.2,0.4,0.6,0.8,1), 
#'               coefficient_matrix=coefficient_matrix)
T_3<-function(z_vector,
              cov_mat,
              cutoff_value,
              coefficient_matrix){
    pvalue<-rep(0,length(cutoff_value))
    for (i in 1:(length(cutoff_value)-1)) {
        stat_TWT<-T_eta(z_vector,cov_mat,cutoff_value[i])
        pvalue[i]<-stats::pgamma(
            (stat_TWT-coefficient_matrix[2,i])/(coefficient_matrix[1,i]),
            shape = coefficient_matrix[3,i]/2,rate = 1/2,lower.tail = FALSE)
    }
    pvalue[length(cutoff_value)]<-chisq_test(z_vector,cov_mat)
    cauchy_statistic<-(1/length(pvalue))*sum(tan((0.5-pvalue)*pi))
    pvalue<-0.5-(atan(cauchy_statistic)/pi)
    return(pvalue)
}