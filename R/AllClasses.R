setClass("anansiWeb",
         slots = c(
           tableY     = "matrix",
           tableX     = "matrix",
           dictionary = "matrix"
           )
         )

setClass("anansiTale",
         slots = c(
           subject   = "character",
           type      = "character",
           estimates = "matrix",
           p.values  = "matrix",
           q.values  = "matrix"
           )
         )
