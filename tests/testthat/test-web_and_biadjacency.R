test_that("kegg wrapper is equivalent", {

  web_g <- weaveWeb(cpd ~ ko, link = list(ec2ko = ec2ko, ec2cpd = ec2cpd))
  web_k <- weaveKEGG(cpd ~ ko)

  expect_identical(web_g, web_k)
})

test_that("Swapping terms in formula is equivalent to transposition", {

  a <- dictionary(weaveWeb(ko ~ cpd, link = kegg_link()))
  b <- dictionary(weaveWeb(cpd ~ ko, link = kegg_link()))

  expect_identical(a, Matrix::t(b))
})
