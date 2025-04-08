test_that("conversion from web to mae to web works", {
    # Make a random anansiWeb
    web1 <- randomWeb()
    # Combine experiments into MultiAssayExperiment object
    mae1 <- asMAE(web1)
    # Back to AnansiWeb
    web2 <- getWeb(x = mae1, tableY = "y", tableX = "x")
    # Once more, my friends
    mae2 <- as(web2, "MultiAssayExperiment")

    expect_identical(web1, web2)
    expect_identical(mae1, mae2)
})
