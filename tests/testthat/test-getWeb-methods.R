test_that("conversion from web to mae to web works", {
    # Make a random anansiWeb
    web <- randomWeb()
    # Combine experiments into MultiAssayExperiment object
    mae <- as(web, "MultiAssayExperiment")
    # Back to AnansiWeb
    outWeb <- getWeb(mae, tableY = "y", tableX = "x")
    expect_identical(web, outWeb)
})
