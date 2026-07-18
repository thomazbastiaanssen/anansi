# Get a listof edges

Get a listof edges

## Usage

``` r
getEdgeList(x, ...)
```

## Arguments

- x:

  input object

- ...:

  additional arguments

## Value

a two-column data.frame that lists the content of each entry in the
input MultiFactor

## Examples

``` r
x <- randomMultiFactor(n_features = 10)
getEdgeList(x)
#>     V1 V2
#> a2b  a  b
#> b2c  b  c
#> c2d  c  d
#> d2e  d  e
#> e2f  e  f
```
