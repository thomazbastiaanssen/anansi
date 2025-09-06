x <- randomWeb()

test_that("AnansiWeb coersion to data.frame and list works", {
    a <- as.data.frame(x)
    b <- do.call(cbind, unname(as.list(x)[c("y", "x", "metadata")]))
    expect_equal(a, b)
})
