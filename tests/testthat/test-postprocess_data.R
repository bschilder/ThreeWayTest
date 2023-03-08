test_that("postprocess_data works", {
  
    dat <- ThreeWayTest::data_matrix_final
    dat2 <- postprocess_data(dat = dat)
    testthat::expect_true(methods::is(dat2,"data.table"))
    testthat::expect_equal(nrow(dat), nrow(dat2))
    testthat::expect_equal(nrow(dat),196036)
    
    
    dat3 <- postprocess_data(dat = dat,
                             agg_var = "GENE")
    testthat::expect_true(methods::is(dat3,"data.table"))
    testthat::expect_gt(nrow(dat), nrow(dat3))
    testthat::expect_equal(nrow(dat3),15829)
})
