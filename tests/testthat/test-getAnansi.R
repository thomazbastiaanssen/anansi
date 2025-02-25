test_that("getAnansi", {
  ### Import libraries ###
  library(mia)
  library(TreeSummarizedExperiment)
  library(MultiAssayExperiment)
  ### Load data ###
  data("FMT_data", package = "anansi")
  ### Prepare objects ###
  metab_se <- SummarizedExperiment(assays = SimpleList(conc = as.matrix(FMT_metab)))
  KO_tse <- TreeSummarizedExperiment(assays = SimpleList(counts = as.matrix(FMT_KOs)))
  keep <- row.names(KO_tse) %in% sort(unique(unlist(anansi_dic)))
  KO_tse <- KO_tse[keep, ]
  KO_tse <- subsetByPrevalent(KO_tse,
    assay.type = "counts",
    prevalence = 0.1
  )
  KO_tse <- transformAssay(KO_tse,
    assay.type = "counts",
    method = "clr",
    pseudocount = TRUE
  )
  coldata <- FMT_metadata
  rownames(coldata) <- coldata$Sample_ID
  coldata <- coldata[match(colnames(KO_tse), rownames(coldata)), ]
  mae <- MultiAssayExperiment(
    experiments = ExperimentList(cpd = metab_se, ko = KO_tse),
    colData = coldata
  )
  ### Check errors and warnings ###
  expect_error(getAnansi(mae),
    "'assay.typeY' must be a valid name of assays(x)",
    fixed = TRUE
  )
  expect_error(getAnansi(mae, experimentY = "wrong_name"),
    "'experimentY' must be numeric or character value specifying experiment in experiment(x)",
    fixed = TRUE
  )
  expect_error(getAnansi(mae,
    experimentY = "cpd",
    experimentX = "ko", assay.typeY = "conc", assay.typeX = "clr",
    return.format = "wrong_input"
  ), class = "error")
  expect_warning(
    getAnansi(mae,
      experimentY = "cpd",
      experimentX = "ko", assay.typeY = "conc", assay.typeX = "clr",
      formula = ~Legend, tableY = 0, tableX = 0
    ),
    "The arguments 'tableY', 'tableX' should not be used, as they are extracted from 'x'",
    fixed = TRUE
  )
  ### Check identity with original anansi output ###
  web <- weaveWeb(
    cpd ~ ko, 
    tableY = t(assay(metab_se, "conc")), tableX = t(assay(KO_tse, "clr")),
    link = kegg_link()
  )
  table1 <- anansi(
    web = web, formula = ~Legend, metadata = FMT_metadata,
    verbose = FALSE
  )
  list1 <- anansi(
    web = web, formula = ~Legend, metadata = FMT_metadata,
    verbose = FALSE, return.format = "list"
  )
  raw1 <- anansi(
    web = web, formula = ~Legend, metadata = FMT_metadata,
    verbose = FALSE, return.format = "raw"
  )

  table2 <- getAnansi(mae,
    experimentY = "cpd",
    experimentX = "ko", assay.typeY = "conc", assay.typeX = "clr",
    formula = ~Legend, translate = TRUE, Y_translation = cpd_translation,
    X_translation = KO_translation, verbose = FALSE
  )
  list2 <- getAnansi(mae,
    experimentY = "cpd",
    experimentX = "ko", assay.typeY = "conc", assay.typeX = "clr",
    formula = ~Legend, translate = TRUE, Y_translation = cpd_translation,
    X_translation = KO_translation, return.format = "list", verbose = FALSE
  )
  raw2 <- getAnansi(mae,
    experimentY = "cpd",
    experimentX = "ko", assay.typeY = "conc", assay.typeX = "clr",
    formula = ~Legend, return.format = "raw", verbose = FALSE
  )

  expect_identical(table1, table2)
  expect_identical(list1, list2)
  expect_identical(raw1, raw2)
})
