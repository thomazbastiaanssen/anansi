test_that("full and disjointed parameters correspond to stats::lm()", {

  tX      <- `colnames<-`(replicate(5, c(scale(rnorm(36)))), letters[1:5])
  tY      <- `colnames<-`(replicate(3, c(scale(rnorm(36)))), LETTERS[1:3])

  d   <- matrix(TRUE, nrow = NCOL(tY), ncol = NCOL(tX), dimnames = list(colnames(tY), colnames(tX)))

  m    <- data.frame(categorical = scale(rnorm(36)),
                            continuous = sample(c("X", "Y", "Z"), size = 36, replace = TRUE))

  web <- new("anansiWeb",
             tableY = tY,
             tableX = tX,
             dictionary = d)


  anansi.res <- anansi:::anansiDiffCor(web, m, formula = ~ categorical, verbose = FALSE)

  a.full.F   <- tell_F(anansi.res$modelfit)[[1]]
  a.disj.F   <- tell_F(anansi.res$disjointed)[[1]]

  a.full.Rsq <- tell_e(anansi.res$modelfit)[[1]]
  a.disj.Rsq <- tell_e(anansi.res$disjointed)[[1]]

  a.full.P   <- tell_P(anansi.res$modelfit)[[1]]
  a.disj.P   <- tell_P(anansi.res$disjointed)[[1]]

  t.full     <- matrix(nrow = NCOL(tY),
                       ncol = NCOL(tX),
                       dimnames = list(colnames(tY),
                                       colnames(tX)))
  t.disj.F <- t.full.F <- t.full.P <- t.disj.P <- t.disj <- t.full


  for(y in 1:NCOL(tY)){
    m$y.val = y.val = tY[,y]

    for(x in 1:NCOL(tX)){
      m$x.val = x.val = tX[,x]

      f.1 <- lm(y.val ~         categorical, data = m)
      f.2 <- lm(y.val ~ x.val * categorical, data = m)
      f.x <- lm(y.val ~ x.val + categorical, data = m)

      df1.f <- f.2$rank - f.1$rank
      df1.x <- f.2$rank - f.x$rank
      df2   <- NROW(f.2$model) - f.2$rank

      ss.1 = sum(residuals(f.1)^2)
      ss.2 = sum(residuals(f.2)^2)
      ss.x = sum(residuals(f.x)^2)

      t.full  [y,x]  <- 1 - (ss.2/ss.1)
      t.full.F[y,x]  <- t.full[y,x] * (df2/df1.f)
      t.full.P[y,x]  <- pf(t.full.F[y,x], df1.f, df2, lower.tail = FALSE)

      t.disj  [y,x]  <- 1 - (ss.2/ss.x)
      t.disj.F[y,x]  <- t.disj[y,x] * (df2/df1.x)
      t.disj.P[y,x]  <- pf(t.disj.F[y,x], df1.x, df2, lower.tail = FALSE)
    }
  }

  #following should be the same, by group
  expect_equal(t.full, a.full.Rsq)
  expect_equal(t.disj, a.disj.Rsq)

  expect_equal(t.full.P, a.full.P)
  expect_equal(t.disj.P, a.disj.P)

  expect_equal(t.full.F, a.full.F)
  expect_equal(t.disj.F, a.disj.F)
})

#
# #Why does ANOVA not give the same F and corresponding p? What is the difference between the two methods?
# #Are there certain conditions that make one method more reasonable to take?
# #
#
# #We can recover same p-value from anova output:
# fval.full <- (1- anova(f.1, f.2)$RSS[2]/anova(f.1, f.2)$RSS[1]) * (anova(f.1, f.2)[2, 1]/anova(f.1, f.2)$Df[2])
# df2.full <- anova(f.1, f.2)[2, 1]
# df1.full <- anova(f.1, f.2)$Df[2]
#
#
# pf(fval.full, df1.full, df2.full, lower.tail = FALSE)
#
# t.full.P[3,5]
# anova(f.1, f.2, scale = 1.0269) #check out scale argument!
# t.full.F[3,5]
#
# t.full.F[3,5]/ anova(f.1, f.2, scale = 1)$F[2]
#
# stats:::anova.lmlist
# ?stat.anova
