test_that("gen_autoregressive works", {
    
    m <- 6
    cov_mat <- gen_autoregressive(m = m,
                                  rho = 0.5)
    testthat::expect_equal(dim(cov_mat),
                           rep(m,2))
})
