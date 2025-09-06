x <- randomWeb()

test_that("getFeaturePairs and mapply work", {
    expect_equal(
        lapply(
            pairwiseApply(
                x,
                FUN = function(x, y) cbind(y, x),
                USE.NAMES = FALSE,
                SIMPLIFY = FALSE
            ),
            unname
        ),
        lapply(getFeaturePairs(x), unname)
    )
})

test_that("show works", {
    expect_null(show(x))
})
