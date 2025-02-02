test_that("getAnansi", {
  ### Import libraries ###
  library(mia)
  library(TreeSummarizedExperiment)
  library(MultiAssayExperiment)
  ### Prepare objects ###
  metab_se <- SummarizedExperiment(assays = SimpleList(conc = as.matrix(FMT_metab)))
  KO_tse <- TreeSummarizedExperiment(assays = SimpleList(counts = as.matrix(FMT_KOs)))
  keep <- row.names(KO_tse) %in% sort(unique(unlist(anansi_dic)))
  KO_tse <- KO_tse[keep, ]
  KO_tse <- subsetByPrevalent(KO_tse,
                              assay.type = "counts",
                              prevalence = 0.1)
  KO_tse <- transformAssay(KO_tse,
                           assay.type = "counts",
                           method = "clr",
                           pseudocount = TRUE)
  coldata <- FMT_metadata
  rownames(coldata) <- coldata$Sample_ID
  coldata <- coldata[match(colnames(KO_tse), rownames(coldata)), ]
  mae <- MultiAssayExperiment(
    experiments = ExperimentList(metabolites = metab_se, functions = KO_tse),
    colData = coldata
  )
  ### Check errors and warnings ###
  expect_error(getAnansi(mae),
      "'assay.type1' must be a valid name of assays(x)", fixed = TRUE)
  expect_error(getAnansi(mae, experiment1 = "wrong_name"),
      "'experiment1' must be numeric or character value specifying experiment in experiment(x)",
      fixed = TRUE)
  expect_error(getAnansi(mae, experiment1 = "metabolites",
      experiment2 = "functions", assay.type1 = "conc", assay.type2 = "clr",
      return.long = "wrong_input"), "'return.long' must be TRUE or FALSE",
      fixed = TRUE)
  expect_warning(getAnansi(mae, experiment1 = "metabolites",
      experiment2 = "functions", assay.type1 = "conc", assay.type2 = "clr",
      formula = ~ Legend, tableY = 0, tableX = 0),
      "The arguments 'tableY', 'tableX' should not be used, as they are extracted from 'x'",
      fixed = TRUE)
  ### Check identity with original anansi output ###
  web <- weaveWebFromTables(tableY = t(assay(metab_se, "conc")), verbose = FALSE,
      tableX = t(assay(KO_tse, "clr")), dictionary = anansi_dic)
  out1 <- anansi(web = web, formula  = ~ Legend, metadata = FMT_metadata, verbose = FALSE)
  long1 <- spinToLong(anansi_output = out1, translate = TRUE, 
      Y_translation = cpd_translation, X_translation = KO_translation)
  out2 <- getAnansi(mae, experiment1 = "metabolites", experiment2 = "functions",
      assay.type1 = "conc", assay.type2 = "clr", formula = ~ Legend,
      return.long = FALSE, verbose = FALSE)
  long2 <- getAnansi(mae, experiment1 = "metabolites", experiment2 = "functions",
      assay.type1 = "conc", assay.type2 = "clr", formula = ~ Legend,
      translate = TRUE, Y_translation = cpd_translation,
      X_translation = KO_translation, verbose = FALSE)
  expect_identical(out1, out2)
  expect_identical(long1, long2)
})
