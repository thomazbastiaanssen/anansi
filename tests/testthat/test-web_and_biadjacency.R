test_that("kegg wrapper is equivalent", {

  generic      <- weaveWeb(cpd ~ ko, link = list(ec2ko, ec2cpd))
  kegg_wrapper <- weaveKEGG(cpd ~ ko)

  expect_identical(generic, kegg_wrapper)
})

test_that("Swapping terms in formula is equivalent to transposition", {

  a <- weaveWeb(ko ~ cpd, link = list(ec2ko, ec2cpd))$dictionary
  b <- weaveWeb(cpd ~ ko, link = list(ec2ko, ec2cpd))$dictionary

  expect_identical(a, Matrix::t(b))
})
