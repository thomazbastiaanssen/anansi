# Apply a function on each pair of features

Apply a function on each pair of features

## Usage

``` r
pairwiseApply(X, ...)
```

## Arguments

- X:

  input object

- ...:

  additional arguments

## Value

a list containing the output of applying the function to each feature
pair. See `?base::mapply()`

## Examples

``` r
web <- randomWeb(10)

# For each feature pair, was the value for x higher than the value for y?
pairwise_gt <- pairwiseApply(
    X = web,
    FUN = function(x, y) x > y,
    MoreArgs = NULL, SIMPLIFY = FALSE, USE.NAMES = TRUE
)

head(pairwise_gt)
#> $x_1y_1
#>  [1] FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE
#> 
#> $x_1y_3
#>  [1]  TRUE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE
#> 
#> $x_1y_5
#>  [1] FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE
#> 
#> $x_1y_6
#>  [1] FALSE FALSE  TRUE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE
#> 
#> $x_1y_8
#>  [1]  TRUE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE
#> 
#> $x_1y_9
#>  [1] FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE
#> 

# Run cor.test() on each pair of features
pairwise_cor <- pairwiseApply(
    X = web,
    FUN = function(x, y) cor.test(x, y),
    MoreArgs = NULL, SIMPLIFY = FALSE, USE.NAMES = TRUE
)

pairwise_cor[1]
#> $x_1y_1
#> 
#>  Pearson's product-moment correlation
#> 
#> data:  x and y
#> t = 0.62311, df = 8, p-value = 0.5506
#> alternative hypothesis: true correlation is not equal to 0
#> 95 percent confidence interval:
#>  -0.4794249  0.7439896
#> sample estimates:
#>       cor 
#> 0.2151447 
#> 
#> 
```
