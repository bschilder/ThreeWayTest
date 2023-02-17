#' quadratic test
#' 
#' This is the function for method used in real data analysis 
#' and simulation result.
#' This is a traditional test, for more information see zhu_Shom
#' input value description: 
#' z_vector is the vector of z_scores
#' estimated_cov is the estimated covariance matrix of these z_scores
#' output value description: 
#' p_value of quadratic test
#' 
#' @export
#' @importFrom MASS ginv
quadratic_test<-function(z_vector, 
                         estimated_cov){
  chi_stat<-t(z_vector)%*%MASS::ginv(estimated_cov)%*%z_vector
  p_value<-stats::pchisq(chi_stat,length(z_vector),lower.tail = FALSE)
  p_value<-as.numeric(p_value)
  return(p_value)
}



#' true_cov_mat_generation_autoregressive
#' 
#' Function for covariance matrix generation.
#' @export
true_cov_mat_generation_autoregressive<-function(m,
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
#' MGAS: this is the method proposed by Van, for more information, please see 
#' Van_MGAS#
#' input value description: 
#' z_vector is the vector of z_scores
#' est_genetic_cov is the estimated correlation matrix of genetic variable
#' est_genetic_cov is the estimated correlation matrix of genetic variable
#' output value description: 
#' p_value of MGAS
MGAS<-function(z_vector,
               est_genetic_cor,
               est_pheno_cor){
  X <- methods::kronecker(est_genetic_cor,est_pheno_cor,FUN = "*")
  beta = c(0.3867,0.0021,-0.1347,-0.0104,0.7276,0.0068)
  p_value<-stats::pchisq(z_vector^2,1,lower.tail = F) # original p-value
  tmp=sort(p_value,index.return=T)
  pj=tmp$x # sorted p-values
  iorder=tmp$ix # index
  length = ncol(est_genetic_cor)*ncol(est_pheno_cor)
  r2=matrix(0,length,length) 
  r2=X[iorder,iorder]
  ro=diag(length)
  for (i1 in 1:length)  {  
    for (i2 in 1:i1) {
      if (i1>i2) {
        er=r2[i1,i2]
        ro[i1,i2]=ro[i2,i1]= beta[1]*er^6+beta[2]*er^5+beta[3]*er^4+beta[4]*er^3+beta[5]*er^2+beta[6]*er
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



#' expM
#' 
#' metaCCA: this is the method proposed by Cichonska, for more information please see cic_metaCCA#
#' code website: https://github.com/aalto-ics-kepaco
#' input value description: 
#' please put in the description
expM <- function(X, 
                 e) {
  
  v	=  La.svd(X)  
  res	=  v$u %*% diag(v$d^e) %*% v$vt
  
  return(res)
} 

#' myCCA
#' 
#' (Add documentation here)
#' @export
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
  p_val = 	pchisq(chi, df, lower.tail = FALSE)
  return(p_val)
  #return( list(r, p_val, a[, 1], b[, 1]) )
}


#' shrinkPSD
#' 
#' (Add documentation here)
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

metaCCA <- function(est_genetic_cov,
                    est_pheno_cov,
                    est_XY_cov,
                    N){
  M<-rbind(cbind(est_genetic_cov,est_XY_cov),cbind(t(est_XY_cov),est_pheno_cov)) #??????
  M_shrink <- shrinkPSD(M,ncol(est_genetic_cov))
  p_metacca<-myCCA(M_shrink[[1]],M_shrink[[2]],M_shrink[[3]],N)
  return(p_metacca[1])
}


#' standard_chisq
#' 
#' This is the function of T_3, sever as a part of new method
#' @export
#' @importFrom MASS ginv
standard_chisq<-function(z_vector,
                         correlation_matrix){
  chi_stat<-t(z_vector)%*%MASS::ginv(correlation_matrix)%*%z_vector
  pvalue<-stats::pchisq(chi_stat,length(z_vector),lower.tail = FALSE)
  pvalue<-as.numeric(pvalue)
  return(pvalue)
}

#' calculate_stat_T3
#' 
#' (Add documentation here)
#' @export
#' @importFrom MASS ginv
calculate_stat_T3<-function(z_vector,
                            correlation_matrix,
                            c_fixed){
  index_c<-ceiling(length(z_vector)*c_fixed)
  z_vector_order<-order(abs(z_vector),decreasing = TRUE)
  ordered_z_vector<-z_vector[z_vector_order]
  ordered_correlation_matrix<-correlation_matrix[z_vector_order,z_vector_order]
  select_vector<-ordered_z_vector[1:index_c]
  select_correlation_matrix<-ordered_correlation_matrix[(1:index_c),(1:index_c)]
  stat<-t(select_vector)%*%MASS::ginv(select_correlation_matrix)%*%select_vector
  stat<-as.numeric(stat)
  return(stat)
}

#' generate_null_distribution_T3
#' 
#' This function generate the null distribution of T_eta, should be run before calculating TWT
#' output generated null distribution for function approximate_distribution_coefficient_estimate_T3
#' Please see Setp 2 in our manuscript for more information
#' m stantds for the dimension, n is the number of sample size
generate_null_distribution_T3<-function(m,
                                        n,
                                        correlation_matrix,
                                        cutoff_value){
  true_mean<-rep(0,m)
  bootstrap_data<-MASS::mvrnorm(n,true_mean,correlation_matrix)
  null_distribution_matrix<-matrix(0,nrow = n, ncol = length(cutoff_value)-1)
  for (i in 1:(length(cutoff_value)-1)) {
    for (j in 1:n) {
      null_distribution_matrix[j,i]<-calculate_stat_T3(bootstrap_data[j,],correlation_matrix,cutoff_value[i])
    }
  }
  return(null_distribution_matrix)
}

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
#' This function input null distribution for generate_null_distribution_T3 and
#' output estimated coefficient. Should be run before calculating TWT
#' Please see Setp 4 in our mauscript for more information 
approximate_distribution_coefficient_estimate_T3<-function(null_distribution_matrix){
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
#' (Add documentation here)
T_3<-function(z_vector,
              correlation_matrix,
              cutoff_value,
              coefficient_matrix){
  pvalue<-rep(0,length(cutoff_value))
  for (i in 1:(length(cutoff_value)-1)) {
    stat_TWT<-calculate_stat_T3(z_vector,correlation_matrix,cutoff_value[i])
    pvalue[i]<-pgamma((stat_TWT-coefficient_matrix[2,i])/(coefficient_matrix[1,i]),shape = coefficient_matrix[3,i]/2,rate = 1/2,lower.tail = F)
  }
  pvalue[length(cutoff_value)]<-standard_chisq(z_vector,correlation_matrix)
  cauchy_statistic<-(1/length(pvalue))*sum(tan((0.5-pvalue)*pi))
  final_pvalue<-0.5-(atan(cauchy_statistic)/pi)
  return(final_pvalue)
}


#' new_method
#' 
#' This is the method of TWT
#' input value description:
#' z_mat is the matrix of z_scores with row stands for number of variants and coloumns standsfor the number of phenotypes.
#' est_genetic_cor stands for the estimated correlation matrix of genetic variants
#' est_pheno_cor stands for the estimated correlation matrix of phenotypes
#' cutoff_value is the set Omega in our paper and the coefficient matrix is calculated based on function approximate_distribution_coefficient_estimate_T3
#' (number of row stands for the number of SNPs)
#' (number of coloumn  stands for the number of Phenotypes)
#' output value description:
#' final p_value and list of final combined 3 pvalues
new_method<-function(z_mat, 
                     est_genetic_cor, 
                     est_pheno_cor, 
                     cutoff_value,
                     coefficient_matrix){
  z_mat<-t(z_mat)
  number_of_snp<-ncol(z_mat)
  number_of_pheno<-nrow(z_mat)
  pvalue_genetic<-rep(NA,number_of_snp)
  pvalue_pheno<-rep(NA,number_of_pheno)
  for (i in 1:number_of_snp) {
    pvalue_genetic[i]<-standard_chisq(z_mat[,i],est_pheno_cor)
  }
  cauchy_statistic_1<-(1/number_of_snp)*sum(tan((0.5-pvalue_genetic)*pi))
  p_pleiotropic<-0.5-(atan(cauchy_statistic_1)/pi)
  for (i in 1:number_of_pheno) {
    pvalue_pheno[i]<-standard_chisq(z_mat[i,],est_genetic_cor)
  }
  cauchy_statistic_2<-(1/number_of_pheno)*sum(tan((0.5-pvalue_pheno)*pi))
  p_genetic_structure<-0.5-(atan(cauchy_statistic_2)/pi)
  z_vec<-as.vector(z_mat)
  est_total_cov_mat<-methods::kronecker(est_genetic_cor,est_pheno_cor)
  p_3<-T_3(z_vec,est_total_cov_mat,cutoff_value,coefficient_matrix)
  final_p_vec<-c(p_1,p_2,p_3)
  cauchy_statistic_final<-(1/3)*sum(tan((0.5-final_p_vec)*pi))
  p_value_final<-0.5-(atan(cauchy_statistic_final)/pi)
  return(list(p_value_final=p_value_final, 
              p_3=p_3,
              p_1=p_genetic_structure,
              p_2=p_pleiotropic))
}
