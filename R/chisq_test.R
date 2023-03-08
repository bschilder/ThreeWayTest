#' chisq_test
#' 
#' This is traditional chisq_test.
#' @param z_vector Column vectorized data matrix with rows represent
#'  phenotype and columns represent genotype.
#' @param cov_mat Estimated covariance matrix of z_vector.
#' @returns A numeric value represent the p value of chi-square test.
#' @importFrom MASS ginv
#' @export
#' @examples
#' z_vector<-MASS::mvrnorm(n = 1,
#'                         mu=rep(0,9),
#'                         Sigma = diag(nrow = 9, ncol = 9))
#' ThreeWayTest::chisq_test(z_vector=z_vector, 
#'                          cov_mat=diag(nrow = 9, ncol = 9))
chisq_test<-function(z_vector,
                     cov_mat){
    chi_stat<-t(z_vector)%*%MASS::ginv(cov_mat)%*%z_vector
    pvalue<-stats::pchisq(chi_stat,length(z_vector),lower.tail = FALSE)
    pvalue<-as.numeric(pvalue)
    return(pvalue)
}