#' gen_autoregressive
#'
#' This is the function for generating autoregressive corvaiance matrix.
#' @param m Dimension of the matrix.
#' @param rho Correlation coefficient.
#' @returns A autoregressive covariance matrix.
#' @export
#' @examples
#' ThreeWayTest::gen_autoregressive(6,0.5)
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


#' MGAS
#'
#' This is the method proposed by Van_MGAS
#' @param z_vector Column vectorized data matrix with rows represent
#'  phenotype and columns represent genotype.
#' @param est_genetic_cor Estimated phenotype correlation matrix.
#' @param est_pheno_cor Estimated genotype correlation matrix.
#' @returns A numeric value represents the p value of MGAS.
#' 
#' @export
#' @importFrom stats pchisq
#' @examples
#' z_vector<-MASS::mvrnorm(1,mu=rep(0,9),Sigma = diag(nrow = 9, ncol = 9))
#' genotype_covariance<-diag(nrow = 3,ncol = 3)
#' phenotype_covariance<-diag(nrow = 3,ncol = 3)
#' ThreeWayTest::MGAS(z_vector=z_vector,
#'                    est_genetic_cor = genotype_covariance, 
#'                    est_pheno_cor =phenotype_covariance)
MGAS<-function(z_vector,est_genetic_cor,est_pheno_cor){
  X <- kronecker(est_genetic_cor,est_pheno_cor,FUN = "*")
  beta = c(0.3867,0.0021,-0.1347,-0.0104,0.7276,0.0068)
  p_value<-stats::pchisq(z_vector^2,1,lower.tail = FALSE) 
  tmp=sort(p_value,index.return=TRUE) 
  pj=tmp$x 
  iorder=tmp$ix
  length = ncol(est_genetic_cor)*ncol(est_pheno_cor)
  r2=matrix(0,length,length) 
  r2=X[iorder,iorder]
  ro=diag(length)
  for (i1 in 1:length)  {  
    for (i2 in 1:i1) {
      if (i1>i2) {
        er=r2[i1,i2]
        ro[i1,i2]=ro[i2,i1]= beta[1]*er^6+beta[2]*er^5+beta[3]*er^4+beta[4]*
            er^3+beta[5]*er^2+beta[6]*er
      }}}
  alllam=eigen(ro[1:length,1:length])$values
  qepj=length 
  for (i1 in 1:length) {
    qepj=qepj-(alllam[i1]>1)*(alllam[i1]-1) }
  
  qej=matrix(c(seq(1,length,1)),length,1,byrow=T)
  for (j in 1:length) { 
    sellam=eigen(ro[1:j,1:j])$values
    id=j
    for (i1 in 1:id) {
      qej[j,1]=qej[j,1]-(sellam[i1]>1)*(sellam[i1]-1)
    }
  }
  pg=matrix(0,length,1)
  for (i in 1:length) {
    pg[i,1]=(qepj/qej[i,1])*pj[i]
  }
  pg=pg[iorder]
  p_mgas<-min(pg)
  return(p_mgas)
}


expM <- function(X, 
                 e) {
  
  v	=  La.svd(X)  
  res	=  v$u %*% diag(v$d^e) %*% v$vt
  
  return(res)
} 

myCCA <- function(C_XX, 
                  C_YY, 
                  C_XY, 
                  N) {
  C_XX = as.matrix(C_XX);  C_YY = as.matrix(C_YY);
  C_XY = as.matrix(C_XY)
  expM_Cxx = expM(C_XX, -0.5)
  expM_Cyy = expM(C_YY, -0.5)
  K    =  expM_Cxx  %*% C_XY  %*%  expM_Cyy 
  # Singular Value Decomposition (SVD)
  svd_out =  svd(K)
  U 		=  svd_out$u #left
  S		=  as.matrix(svd_out$d)
  V		=  svd_out$v
  # S is a vector, create a diagonal matrix
  temp = diag(length(S))
  diag(temp) = S
  
  min_r 	=  min( dim(C_XX)[1], dim(C_YY)[1] )
  a 		=  array( dim = c( dim(C_XX)[1], min_r) )  	# canonical weights
  b 		=  array( dim = c( dim(C_YY)[1], min_r) )
  r 		=  array( dim = c(1, min_r) )             	# canonical correlation
  wilks 	=  array( dim = c(1, min_r) )         		# Wilk's lambda
  lambdas	=  array( dim = c(length(S), length(S)) )
  lambdas =  temp^2
  df 		=  array( dim = c(1, min_r) )             	# degrees of freedom
  p 		=  dim(C_XX)[1]
  q 		=  dim(C_YY)[1]
  k		=  0;
  for ( i in 1:min_r ) {
    a[,i]	=  expM_Cxx  %*%  U[,i]
    b[,i] 	=  expM_Cyy  %*%  V[,i]
    r[1,i]	=  ( t(a[,i])%*%C_XY%*%b[,i] ) /
      ( sqrt( t(a[,i])%*%C_XX%*%a[,i] )  *
          sqrt( t(b[,i])%*%C_YY%*%b[,i] ) )
    # Wilk's lambda
    prodd = 1
    for (j in i:min_r) {
      prod_temp 	=  1 - lambdas[j,j]
      prodd 	=  prodd * prod_temp
    }
    wilks[i] 	=  prodd
    df[i] 		=  (p-k)*(q-k)       # df for Bartlett Chi-square
    k 			=  k+1
  }
  # Bartlett Chi-square approximation to Wilk's Lambda
  chi 	= 	-( (N-1) - 0.5*(p+q+1) )  *  log(wilks)
  # H0: all canonical correlations = 0
  # p-value
  p_val = 	stats::pchisq(chi, df, lower.tail = FALSE)
  return(p_val)
  #return( list(r, p_val, a[, 1], b[, 1]) )
}


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

#' coefficient_estimate
#' 
#' This function estimate the coeffcient of T3 with null distribution 
#' generated by generate_null_distribution_T3.
#' Please see Setp 2 in our manuscript for more information.
#' @param null_distribution Generated by generate_null_distribution_T3
#' @returns List of coefficient of alpha, beta and d.
#' @export
#' @examples
#' null_distribution<-ThreeWayTest::generate_null_distribution_T3(m=6,n=1000,
#' cov_mat=diag(nrow = 6, ncol = 6), cutoff_value=c(0.2,0.4,0.6,0.8,1))
#' ThreeWayTest::coefficient_estimate(null_distribution[,1])
coefficient_estimate<-function(null_distribution){
  n<-length(null_distribution)
  K1<-(1/n)*sum(null_distribution)
  K2<-(1/n)*sum((null_distribution-K1)^2)
  K3<-(1/n)*sum((null_distribution-K1)^3)
  K4<-(1/n)*sum((null_distribution-K1)^4)-3*(K2^2)
  alpha<-K3/(4*K2)
  beta<-K1-2*(K2^2)/K3
  d<-8*(K2^3)/(K3^2)
  return(list(alpha=alpha,beta=beta,d=d))
}

#' approximate_distribution_coefficient_estimate_T3
#' 
#' This function estimate the coeffcient of T3 with null distribution 
#' generated by  generate_null_distribution_T3.
#' Please see Setp 2 in our manuscript for more information.
#' @param null_distribution_matrix Generated by generate_null_distribution_T3
#' @returns Coefficient matrix of alpha, beta and d.
#' @export
#' @examples 
#' null_distribution<-ThreeWayTest::generate_null_distribution_T3(m=6,n=1000,
#' cov_mat=diag(nrow = 6, ncol = 6), cutoff_value=c(0.2,0.4,0.6,0.8,1))
#' coefficient_matrix <- 
#'     ThreeWayTest::approximate_distribution_coefficient_estimate_T3(
#'     null_distribution_matrix = null_distribution)
approximate_distribution_coefficient_estimate_T3<-
    function(null_distribution_matrix){
  coefficient_matrix<-matrix(0,nrow = 3,ncol = ncol(null_distribution_matrix))
  for (i in 1:ncol(null_distribution_matrix)) {
    coefficient<-coefficient_estimate(null_distribution_matrix[,i])
    coefficient_matrix[1,i]<-coefficient$alpha
    coefficient_matrix[2,i]<-coefficient$beta
    coefficient_matrix[3,i]<-coefficient$d
  }
  return(coefficient_matrix)
}

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
#' z_vector<-MASS::mvrnorm(1,mu=rep(0,10), Sigma = diag(nrow = 10, ncol = 10))
#' null_distribution<-ThreeWayTest::generate_null_distribution_T3(m=10,n=1000,
#' cov_mat=diag(nrow = 10,ncol= 10), cutoff_value=c(0.2,0.4,0.6,0.8,1))
#' coefficient_matrix<-
#'     ThreeWayTest::approximate_distribution_coefficient_estimate_T3(
#'         null_distribution_matrix = null_distribution)
#' ThreeWayTest::T_3(z_vector=z_vector, cov_mat=diag(nrow = 10, ncol = 10), 
#' cutoff_value=c(0.2,0.4,0.6,0.8,1), coefficient_matrix=coefficient_matrix)
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
#' z_mat<-MASS::mvrnorm(5,mu=rep(0,5), Sigma = diag(nrow = 5, ncol = 5))
#' null_distribution<-ThreeWayTest::generate_null_distribution_T3(m=25,
#' n=1000,cov_mat=diag(nrow = 25,ncol= 25), cutoff_value=c(0.2,0.4,0.6,0.8,1))
#' coefficient_matrix<-
#'     ThreeWayTest::approximate_distribution_coefficient_estimate_T3(
#'         null_distribution_matrix = null_distribution)
#' ThreeWayTest::TWT(z_mat=z_mat, est_genetic_cor=diag(nrow = 5, ncol = 5), 
#' est_pheno_cor=diag(nrow = 5, ncol = 5), cutoff_value=c(0.2,0.4,0.6,0.8,1), 
#' coefficient_matrix=coefficient_matrix)
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
