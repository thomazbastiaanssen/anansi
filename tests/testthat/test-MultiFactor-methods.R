x <- MultiFactor(kegg_link(), drop.unmatched = FALSE)

test_that("Delayed shedding levels works; MultiFactor returns MultiFactor", {
    expect_identical(
        MultiFactor(x),
        MultiFactor(x)
    )
})

test_that("Dropping levels works", {
    x <- randomMultiFactor(n_features = 10)

    expect_error(
        droplevels(
            x,
            select = list(a = "a_001"),
            exclude = list(d = c("d_001"))
        ),
        regexp = "Only one of 'exclude' and 'select' may be provided"
    )
    y <- droplevels(x, select = list(a = "a_010"))
    z <- droplevels(x, exclude = list(a = paste0("a_00", seq(1, 9))))

    expect_identical(y, z)

    expect_identical(dim(y), c(5L, 6L))
    expect_identical(names(names(y)), rownames(y))
    expect_identical(y, y[])
})

x <- randomMultiFactor(n_features = 10, n_types = 3)

test_that("MultiFactor indexing works", {
    expect_identical(
        x,
        asMultiFactor(x[["b"]], levels = levels(x))
    )

    expect_identical(
        x,
        MultiFactor(c(x[["a"]], x[["c"]]), levels = levels(x))
    )
})

test_that("MultiFactor get/set works", {
    # two-way equivalence
    expect_identical(x[["a"]], x[["a"]] <- x[, c("a", "b")][1])
    expect_identical(x[, c("a", "b")][1], x[, c("a", "b")][1] <- x[["a"]])

    expect_identical(x@map, x@map <- x@map)
})


test_that("show works", {
    expect_null(show(x))
})
