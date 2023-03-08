#' TWT
#'
#' This is the method of TWT.
#' @param z_mat is the matrix of z_scores with row number of rows stands for
#' the number of phenotypes and number of columns
#'  stands for the number of variants.
#' @param est_genetic_cor Estimated correlation matrix of genetic variants.
#' @param est_pheno_cor Estimated correlation matrix of phenotypes.
#' @param cutoff_value Set Omega.
#' @param coefficient_matrix Calculated based on function 
#' \link[ThreeWayTest]{approximate_distribution_coefficient_estimate_T3}.
#' @returns p_value_final P_value of TWT.
#' @returns p_1 P_value of T_1.
#' @returns p_2 P_value of T_2.
#' @returns p_3 P_value of T_3.
#' 
#' @export
#' @examples
#' z_mat<-MASS::mvrnorm(n = 5,
#'                      mu=rep(0,5), 
#'                      Sigma = diag(nrow = 5, ncol = 5))
#' null_distribution<-generate_null_distribution_T3(
#'     m=25,
#'     n=1000,
#'     cov_mat=diag(nrow = 25,ncol= 25), 
#'     cutoff_value=c(0.2,0.4,0.6,0.8,1))
#' coefficient_matrix<- approximate_distribution_coefficient_estimate_T3(
#'         null_distribution_matrix = null_distribution)
#' res <- TWT(z_mat=z_mat,
#'            est_genetic_cor=diag(nrow = 5, ncol = 5),
#'            est_pheno_cor=diag(nrow = 5, ncol = 5), 
#'            cutoff_value=c(0.2,0.4,0.6,0.8,1),
#'            coefficient_matrix=coefficient_matrix)
TWT<-function(z_mat, 
              est_genetic_cor, 
              est_pheno_cor, 
              cutoff_value,
              coefficient_matrix){
  number_of_snp<-ncol(z_mat)
  number_of_pheno<-nrow(z_mat)
  pvalue_genetic<-apply(z_mat,2,chisq_test, cov_mat=est_pheno_cor)
  cauchy_statistic_1<-(1/number_of_snp)*sum(tan((0.5-pvalue_genetic)*pi))
  p_pleiotropic<-0.5-(atan(cauchy_statistic_1)/pi)
  pvalue_pheno<-apply(z_mat,1,chisq_test, cov_mat=est_genetic_cor)
  cauchy_statistic_2<-(1/number_of_pheno)*sum(tan((0.5-pvalue_pheno)*pi))
  p_genetic_structure<-0.5-(atan(cauchy_statistic_2)/pi)
  z_vec<-as.vector(z_mat)
  est_total_cov_mat<-methods::kronecker(est_genetic_cor,est_pheno_cor)
  p_3<-T_3(z_vec,est_total_cov_mat,cutoff_value,coefficient_matrix)
  final_p_vec<-c(p_genetic_structure,p_pleiotropic,p_3)
  cauchy_statistic_final<-(1/3)*sum(tan((0.5-final_p_vec)*pi))
  p_value_final<-0.5-(atan(cauchy_statistic_final)/pi)
  return(list(p_value_final=p_value_final, 
              p_1=p_genetic_structure,
              p_2=p_pleiotropic, p_3=p_3))
}
