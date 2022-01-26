setClass("anansiWeb",
         slots = c(
           tableY     = "matrix",
           tableX     = "matrix",
           dictionary = "matrix"
           )
         )

setClass("anansiTale",
         slots = c(
           type      = "character",
           estimates = "matrix",
           p.values  = "matrix",
           q.values  = "matrix"
           )
         )

