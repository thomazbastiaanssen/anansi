test_that("full and disjointed parameters correspond to stats::lm()", {
  tX <- `colnames<-`(replicate(5, c(scale(rnorm(36)))), letters[1:5])
  tY <- `colnames<-`(replicate(3, c(scale(rnorm(36)))), LETTERS[1:3])

  d <- matrix(TRUE, nrow = NCOL(tY), ncol = NCOL(tX), dimnames = list(colnames(tY), colnames(tX)))

  m <- data.frame(
    x = 1,
    continuous = scale(rnorm(36)),
    categorical = sample(c("X", "Y", "Z"), size = 36, replace = TRUE)
  )

  web <- new("anansiWeb",
    tableY = tY,
    tableX = tX,
    dictionary = d
  )

  anansi.res <- anansi(web = web, metadata = m,
                       formula = ~categorical, verbose = FALSE)@output@model_results

  a.full.F <- tell_F(anansi.res$modelfit)
  a.disj.F <- tell_F(anansi.res$disjointed)

  a.full.Rsq <- tell_e(anansi.res$modelfit)
  a.disj.Rsq <- tell_e(anansi.res$disjointed)

  a.full.P <- tell_P(anansi.res$modelfit)
  a.disj.P <- tell_P(anansi.res$disjointed)

  t.full <- matrix(
    nrow = NCOL(tY),
    ncol = NCOL(tX),
    dimnames = list(
      colnames(tY),
      colnames(tX)
    )
  )
  t.disj.F <- t.full.F <- t.full.P <- t.disj.P <- t.disj <- t.full


  for (y in 1:NCOL(tY)) {
    m$y.val <- y.val <- tY[, y]

    for (x in 1:NCOL(tX)) {
      m$x.val <- x.val <- tX[, x]

      f.1 <- lm(y.val ~ categorical, data = m)
      f.2 <- lm(y.val ~ x.val * categorical, data = m)
      f.x <- lm(y.val ~ x.val + categorical, data = m)

      df1.f <- f.2$rank - f.1$rank
      df1.x <- f.2$rank - f.x$rank
      df2 <- NROW(f.2$model) - f.2$rank

      ss.1 <- sum(residuals(f.1)^2)
      ss.2 <- sum(residuals(f.2)^2)
      ss.x <- sum(residuals(f.x)^2)

      t.full[y, x] <- 1 - (ss.2 / ss.1)
      t.full.F[y, x] <- oddify(t.full[y, x]) * (df2 / df1.f)
      t.full.P[y, x] <- pf(t.full.F[y, x], df1.f, df2, lower.tail = FALSE)

      t.disj[y, x] <- 1 - (ss.2 / ss.x)
      t.disj.F[y, x] <- oddify(t.disj[y, x]) * (df2 / df1.x)
      t.disj.P[y, x] <- pf(t.disj.F[y, x], df1.x, df2, lower.tail = FALSE)
    }
  }

  # following should be the same, by group
  expect_equal(t.full, a.full.Rsq)
  expect_equal(t.disj, a.disj.Rsq)

  expect_equal(t.full.P, a.full.P)
  expect_equal(t.disj.P, a.disj.P)

  expect_equal(t.full.F, a.full.F)
  expect_equal(t.disj.F, a.disj.F)


  ### And test multiple testing through random intercepts with Error() notation

  anansi.res <- anansi(web = web, metadata = m,
                       formula = ~continuous + Error(categorical),
                       verbose = FALSE)@output@model_results

  a.full.F <- tell_F(anansi.res$modelfit)
  a.disj.F <- tell_F(anansi.res$disjointed)

  a.full.Rsq <- tell_e(anansi.res$modelfit)
  a.disj.Rsq <- tell_e(anansi.res$disjointed)

  a.full.P <- tell_P(anansi.res$modelfit)
  a.disj.P <- tell_P(anansi.res$disjointed)

  t.full <- matrix(
    nrow = NCOL(tY),
    ncol = NCOL(tX),
    dimnames = list(
      colnames(tY),
      colnames(tX)
    )
  )
  t.disj.F <- t.full.F <- t.full.P <- t.disj.P <- t.disj <- t.full

  for (y in 1:NCOL(tY)) {
    m$y.val <- y.val <- tY[, y]

    for (x in 1:NCOL(tX)) {
      m$x.val <- x.val <- tX[, x]

      f.1 <- lm(y.val ~ categorical + continuous, data = m)
      f.2 <- lm(y.val ~ categorical + x.val * continuous, data = m)
      f.x <- lm(y.val ~ categorical + x.val + continuous, data = m)

      df1.f <- f.2$rank - f.1$rank
      df1.x <- f.2$rank - f.x$rank
      df2 <- NROW(f.2$model) - f.2$rank

      ss.1 <- sum(residuals(f.1)^2)
      ss.2 <- sum(residuals(f.2)^2)
      ss.x <- sum(residuals(f.x)^2)

      t.full[y, x] <- 1 - (ss.2 / ss.1)
      t.full.F[y, x] <- oddify(t.full[y, x]) * (df2 / df1.f)
      t.full.P[y, x] <- pf(t.full.F[y, x], df1.f, df2, lower.tail = FALSE)

      t.disj[y, x] <- 1 - (ss.2 / ss.x)
      t.disj.F[y, x] <- oddify(t.disj[y, x]) * (df2 / df1.x)
      t.disj.P[y, x] <- pf(t.disj.F[y, x], df1.x, df2, lower.tail = FALSE)
    }
  }

  # following should be the same, by group
  expect_equal(t.full, a.full.Rsq)
  expect_equal(t.disj, a.disj.Rsq)

  expect_equal(t.full.P, a.full.P)
  expect_equal(t.disj.P, a.disj.P)

  expect_equal(t.full.F, a.full.F)
  expect_equal(t.disj.F, a.disj.F)
})
