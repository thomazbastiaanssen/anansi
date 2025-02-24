test_that("plotAnansi", {
  library(ggplot2)
  data(dictionary)
  data(FMT_data)
  KOs <- floor(FMT_KOs)
  KOs <- apply(KOs, c(1, 2), function(x) as.numeric(as.character(x)))
  KOs <- KOs[apply(KOs == 0, 1, sum) <= (ncol(KOs) * 0.90), ]
  KOs <- KOs[row.names(KOs) %in% sort(unique(unlist(anansi_dic))), ]
  KOs.exp <- clr_c(KOs)
  t1 <- t(FMT_metab)
  t2 <- t(KOs.exp)
  web <- weaveWebFromTables(
    tableY = t1,
    tableX = t2,
    dictionary = anansi_dic
  )
  out <- anansi(
    web = web,
    formula = ~Legend,
    groups = "Legend",
    metadata = FMT_metadata,
    adjust.method = "BH",
  )
  # Check arguments
  expect_no_error(plotAnansi(out))
  expect_error(plotAnansi(out, association.type = "wrong"))
  expect_error(plotAnansi(out, model.var = "wrong"))
  expect_warning(plotAnansi(out, association.type = "full",
      model.var = "Legend"),
      "'model.var' is ignored when 'association type' is set to full")
  expect_warning(plotAnansi(out, signif.threshold = 0.05),
      "'signif.threshold' is ignored when 'association type' is not defined")
  expect_error(plotAnansi(out, association.type = "emergent"),
      "'model.var' must specify a variable of the anansi model when 'association type' is set to emergent",
      fixed = TRUE)
  expect_error(plotAnansi(out, association.type = "disjointed",
      model.var = "Legend", shape_by = "wrong"),
      "'shape_by' must be a character string specifying the name of a 'groups' term used in the original anansi call",
      fixed = TRUE)
  expect_no_error(plotAnansi(out, association.type = "disjointed",
      model.var = "Legend", shape_by = "Legend"))
  # Check output plot
  p <- plotAnansi(out, association.type = "emergent", model.var = "Legend",
      fill_by = "group")
  expect_length(p$guides$guides, 2)
  expect_equal(dim(p$data), c(200, 8))
  expect_false(any(is.na(p$data[["fill"]])))
})
