x <- randomWeb()

test_that("getFeaturePairs and mapply work", {
    expect_equal(
        lapply(
            mapply(x, FUN = function(x, y) cbind(y, x), SIMPLIFY = FALSE),
            unname
        ),
        lapply(getFeaturePairs(x), unname)
    )
})

test_that("show works", {
    expect_null(show(x))
})

test_that("AnansiWeb coersion to data.frame and list works", {
    a <- as.data.frame(x)
    b <- do.call(cbind, unname(as.list(x)[c("y", "x", "metadata")]))
    expect_equal(a, b)
})

test_that("AnansiWeb metadata is robust", {
    expect_identical(
        metadata(x, simplify = TRUE) <- list(metadata(x)),
        metadata(x, simplify = FALSE) <- metadata(x, simplify = FALSE)
    )
})
