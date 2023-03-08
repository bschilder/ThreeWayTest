#' metaCCA
#' 
#' This is the main function of metaCCA.
#' @param est_genetic_cov Estimated genotype correlation matrix.
#' @param est_pheno_cov Estimated phenotype correlation matrix.
#' @param est_XY_cov Estimated correlation matrix between 
#' genotype and phenotype.
#' @param N Sample size.
#' @returns A numeric value represents the p value of metaCCA.
#' @export
#' @examples
#' est_genetic_cov <- diag(nrow = 2,ncol = 2)
#' est_pheno_cov <- diag(nrow = 2,ncol = 2)
#' est_XY_cov <- matrix(rnorm(rep(0,4),diag(nrow = 2, ncol = 2)),2,2)
#' N <- 1000
#' ThreeWayTest::metaCCA(est_genetic_cov,est_pheno_cov,est_XY_cov,N)
metaCCA <- function(est_genetic_cov,
                    est_pheno_cov,
                    est_XY_cov,
                    N){
    M<-rbind(cbind(est_genetic_cov,est_XY_cov),
             cbind(t(est_XY_cov),est_pheno_cov))
    M_shrink <- shrinkPSD(M,ncol(est_genetic_cov))
    p_metacca<-myCCA(M_shrink[[1]],M_shrink[[2]],M_shrink[[3]],N)
    return(p_metacca[1])
}