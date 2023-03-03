test_that("TWT works", {
  z_mat<-MASS::mvrnorm(5,mu=rep(0,5), Sigma = diag(nrow = 5, ncol = 5))
  null_distribution<-ThreeWayTest::generate_null_distribution_T3(m=25,n=1000,
                                                                 cov_mat=diag(nrow = 25,ncol= 25),
                                                                 cutoff_value=c(0.2,0.4,0.6,0.8,1))
  coefficient_matrix<-ThreeWayTest::approximate_distribution_coefficient_estimate_T3(null_distribution)
  result<-ThreeWayTest::TWT(z_mat=z_mat, 
                            est_genetic_cor=diag(nrow = 5, ncol= 5), 
                            est_pheno_cor=diag(nrow = 5, ncol = 5),
                            cutoff_value=c(0.2,0.4,0.6,0.8,1), 
                            coefficient_matrix=coefficient_matrix)
  testthat::expect_lte(result$p_value_final, 1)
  testthat::expect_gte(result$p_value_final, 0)
  
})

test_that("metaCCA works", {
  est_genetic_cov <- diag(nrow = 2,ncol = 2)
  est_pheno_cov <- diag(nrow = 2,ncol = 2)
  est_XY_cov <- matrix(rnorm(rep(0,4),diag(nrow = 2, ncol = 2)),2,2)
  N <- 1000
  result<-ThreeWayTest::metaCCA(est_genetic_cov,est_pheno_cov,est_XY_cov,N)
  testthat::expect_lte(result, 1)
  testthat::expect_gte(result, 0)
})

test_that("MGAS works", {
  z_vector<-MASS::mvrnorm(1,mu=rep(0,9),Sigma = diag(nrow = 9, ncol = 9))
  genotype_covariance<-diag(nrow = 3,ncol = 3)
  phenotype_covariance<-diag(nrow = 3,ncol = 3)
  result<-ThreeWayTest::MGAS(z_vector=z_vector, 
                             est_genetic_cor = genotype_covariance, 
                             est_pheno_cor =phenotype_covariance)
  testthat::expect_lte(result, 1)
  testthat::expect_gte(result, 0)
})
