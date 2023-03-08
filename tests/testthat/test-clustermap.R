test_that("clustermap works", {
  
    run_tests <- function(cm,
                          annot_vars){
        testthat::expect_true(methods::is(cm$plot,"plotly"))
        testthat::expect_true(methods::is(cm$data$X,"matrix"))
        if(is.null(annot_vars)){
            testthat::expect_null(cm$data$row_side_colors)
        } else {
            testthat::expect_true(methods::is(cm$data$row_side_colors,
                                              "data.table"))
        }
    }
    
    #### As raw matrix ####
    dat <- ThreeWayTest::data_matrix_final
    cm1 <- clustermap(dat = dat)
    run_tests(cm=cm1,
              annot_vars=eval(formals(clustermap)$annot_vars))
    testthat::expect_false(nrow(cm1$data$X)==ncol(cm1$data$X))
    
    #### As correlation matrix ####
    annot_vars <- NULL
    cm2 <- clustermap(dat = dat, 
                      as_cor = TRUE,
                      annot_vars = annot_vars)
    run_tests(cm=cm2, 
              annot_vars=annot_vars)
    testthat::expect_equal(nrow(cm2$data$X),
                           ncol(cm2$data$X))
})
