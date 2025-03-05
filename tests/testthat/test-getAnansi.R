test_that("getAnansi", {
  ### Import libraries ###
  library(mia)
  library(TreeSummarizedExperiment)
  library(MultiAssayExperiment)
  ### Load data ###
  ###

  data("FMT_data", package = "anansi")

  # Make a random anansiWeb
  web <- randomWeb()

  # Combine experiments into MultiAssayExperiment object
  mae <- as(web, "MultiAssayExperiment")

  expect_error(getAnansi(mae, tableY = "wrong_name"),
    "'tableY' must be numeric or character value specifying experiment in experiment(x)",
    fixed = TRUE
  )
  expect_error(getAnansi(mae,
    tableY = "y",
    tableX = "x",
    return.format = "wrong_input"
  ), class = "error")
  expect_warning(
    getAnansi(mae,
      tableY = "y",
      tableX = "x",
      formula = ~ cat_ab, web = 0
    ),
    "The arguments 'web' should not be used, as they are extracted from 'x'",
    fixed = TRUE
  )
  ### Check identity with original anansi output ###
  web <- randomWeb(n_samples = 100)

  table1 <- anansi(
    web = web, formula = ~ cat_XYZ,
    verbose = FALSE
  )
  list1 <- anansi(
    web = web, formula = ~cat_XYZ, metadata = FMT_metadata,
    verbose = FALSE, return.format = "list"
  )
  raw1 <- anansi(
    web = web, formula = ~cat_XYZ, metadata = FMT_metadata,
    verbose = FALSE, return.format = "raw"
  )

  table2 <- getAnansi(mae,
                      tableY = "y",
                      tableX = "x",
                      formula = ~cat_XYZ, verbose = FALSE)
  list2 <- getAnansi(mae,
                     tableY = "y",
                     tableX = "x",
                     formula = ~cat_XYZ,
                     return.format = "list", verbose = FALSE
  )
  raw2 <- getAnansi(mae,
                    tableY = "y",
                    tableX = "x",
                    formula = ~cat_XYZ,
                    return.format = "raw", verbose = FALSE
  )

  expect_identical(table1, table2)
  expect_identical(list1, list2)
  expect_identical(raw1, raw2)
})
