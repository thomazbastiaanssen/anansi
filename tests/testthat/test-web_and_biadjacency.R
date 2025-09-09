test_that("kegg wrapper is equivalent", {
    data("ec2cpd", package = "anansi", envir = environment())
    data("ec2ko", package = "anansi", envir = environment())
    expect_warning(
        web_g <- weaveWeb(cpd ~ ko, link = list(ec2ko = ec2ko, ec2cpd = ec2cpd)),
        "Argument `metadata` not provided; Please validate sample ID order.")
    expect_warning(
    web_k <- weaveKEGG(cpd ~ ko),
    "Argument `metadata` not provided; Please validate sample ID order.")
    expect_identical(web_g, web_k)
})

test_that("Swapping terms in formula is equivalent to transposition", {
    expect_warning(
        a <- weaveWeb(ko ~ cpd, link = kegg_link())@dictionary,
        "Argument `metadata` not provided; Please validate sample ID order.")
    expect_warning(
    b <- weaveWeb(cpd ~ ko, link = kegg_link())@dictionary,
    "Argument `metadata` not provided; Please validate sample ID order.")

    expect_identical(a, Matrix::t(b))
})
