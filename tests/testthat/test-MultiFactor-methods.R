x <- MultiFactor(kegg_link())

test_that("MultiFactor on a MultiFactor returns itself", {
    expect_identical(
        x,
        MultiFactor(x)
    )
})


test_that("MultiFactor indexing works", {
    expect_identical(
        x,
        asMultiFactor(x[["ec"]], levels = levels(x))
    )

    expect_identical(
        x,
        MultiFactor(c(x[["ko"]], x[["cpd"]]), levels = levels(x))
    )
})
