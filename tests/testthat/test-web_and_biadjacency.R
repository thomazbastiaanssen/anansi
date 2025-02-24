test_that("kegg wrapper is equivalent", {

  generic      <- web(cpd ~ ko, link = list(ec2ko, ec2cpd))
  kegg_wrapper <- kegg_web(cpd ~ ko)

  expect_identical(generic, kegg_wrapper)
})

test_that("Swapping terms in formula is equivalent to transposition", {

  a <- get_dict(web(ko ~ cpd, link = list(ec2ko, ec2cpd)))
  b <- get_dict(web(cpd ~ ko, link = list(ec2ko, ec2cpd)))

  expect_identical(a, Matrix::t(b))
})
