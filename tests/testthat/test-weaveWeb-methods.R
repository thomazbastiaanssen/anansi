test_that("weaveWeb mae and tse methods work", {
    # Combine experiments into MultiAssayExperiment object
    web <- randomWeb(n_samples = 15, n_reps = 1)
    web@metadata$cat_XYZ <- rep(c("X", "Y", "Z"), 5)
    mae <- asMAE(web)

    web2 <- weaveWeb(x = mae, tableY = "y", tableX = "x")
    # Once more, my friends
    mae2 <- asMAE(web)

    expect_identical(web, web2)
    expect_identical(mae, mae2)

    expect_error(
        weaveWeb(mae, tableY = "wrong_name"),
        "'tableY' must be numeric or character value specifying experiment in experiment(x)",
        fixed = TRUE
    )

    expect_error(
        weaveWeb(
            mae,
            tableY = "y",
            tableX = "x",
            web = 0
        ),
        "unused argument (web = 0)",
        fixed = TRUE
    )
    ### Check identity with original anansi output ###

    table1 <- anansi(
        web = web,
        formula = ~cat_XYZ,
        verbose = FALSE
    )
    list1 <- anansi(
        web = web,
        formula = ~cat_XYZ,
        verbose = FALSE,
        return.format = "list"
    )
    raw1 <- anansi(
        web = web,
        formula = ~cat_XYZ,
        verbose = FALSE,
        return.format = "raw"
    )

    table.mae <- weaveWeb(
        mae,
        tableY = "y",
        tableX = "x"
    ) |> anansi(
        formula = ~cat_XYZ,
        verbose = FALSE
    )
    list.mae <- weaveWeb(
        mae,
        tableY = "y",
        tableX = "x"
    ) |> anansi(
        formula = ~cat_XYZ,
        return.format = "list",
        verbose = FALSE
    )
    raw.mae <- weaveWeb(
        mae,
        tableY = "y",
        tableX = "x"
    ) |> anansi(
        formula = ~cat_XYZ,
        return.format = "raw",
        verbose = FALSE
    )

    expect_identical(table1, table.mae)
    expect_identical(list1, list.mae)
    expect_identical(raw1, raw.mae)

    tse <- asTSE(web)

    table.tse <- weaveWeb(
        tse,
        tableY = "y",
        tableX = "x"
    ) |> anansi(
        formula = ~cat_XYZ,
        verbose = FALSE
    )
    list.tse <- weaveWeb(
        tse,
        tableY = "y",
        tableX = "x"
    ) |> anansi(
        formula = ~cat_XYZ,
        return.format = "list",
        verbose = FALSE
    )
    raw.tse <- weaveWeb(
        tse,
        tableY = "y",
        tableX = "x"
    ) |> anansi(
        formula = ~cat_XYZ,
        return.format = "raw",
        verbose = FALSE
    )


    expect_identical(table1, table.tse)
    expect_identical(list1, list.tse)
    expect_identical(raw1, raw.tse)
})
