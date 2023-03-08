test_that("get_full_genotype works", {
  
    full_genotype <- get_full_genotype()
    testthat::expect_true(methods::is(full_genotype,"data.table"))
    testthat::expect_equal(nrow(full_genotype), 196036)
    testthat::expect_equal(ncol(full_genotype), 2510)
    
    
    testthat::expect_error(
        get_full_genotype(.token = "")
    )
})
