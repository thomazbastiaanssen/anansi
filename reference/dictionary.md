# Get dictionary

Get dictionary

## Usage

``` r
dictionary(x, ...)
```

## Arguments

- x:

  input object

- ...:

  additional arguments

## Value

`dictionary` slot of x

## Examples

``` r
x <- randomWeb(10)
dictionary(x)
#> 12 x 8 sparse Matrix of class "ngCMatrix"
#>       x
#> y      x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8
#>   y_1    |   .   |   |   .   .   |   |
#>   y_2    |   |   .   |   .   .   |   |
#>   y_3    |   |   .   .   .   .   .   .
#>   y_4    |   |   .   |   |   |   .   |
#>   y_5    .   .   |   |   |   |   |   .
#>   y_6    .   |   .   .   |   |   .   |
#>   y_7    .   |   .   .   |   .   .   .
#>   y_8    |   |   .   .   .   |   .   .
#>   y_9    |   |   .   |   .   |   .   .
#>   y_10   .   .   |   .   |   .   |   |
#>   y_11   |   .   |   |   .   .   .   |
#>   y_12   |   .   .   .   |   .   |   |
```
