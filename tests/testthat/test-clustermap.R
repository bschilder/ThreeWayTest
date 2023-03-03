test_that("clustermap works", {
  
    dat <- ThreeWayTest::data_matrix_final
    cm <- clustermap(dat = dat)
    testthat::expect_true(methods::is(cm$plot,"plotly"))
    testthat::expect_true(methods::is(cm$data$X,"matrix"))
    testthat::expect_true(methods::is(cm$data$row_side_colors,"data.table"))
})
