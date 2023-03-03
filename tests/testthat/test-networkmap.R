test_that("networkmap works", {
  
    dat <- ThreeWayTest::data_matrix_final
    nm <- networkmap(dat = dat)
    testthat::expect_true(methods::is(nm$plot,"gg"))
    testthat::expect_true(methods::is(nm$data$graph,"tbl_graph"))
    testthat::expect_true(methods::is(nm$data$nodes,"data.table"))
    testthat::expect_true(methods::is(nm$data$edges,"data.table"))
})
