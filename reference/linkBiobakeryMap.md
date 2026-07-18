# make a link data.frame for biobakery mapping files input

make a link data.frame for biobakery mapping files input

## Usage

``` r
linkBiobakeryMap(map)
```

## Arguments

- map:

  `Character`, result from
  [`readLines()`](https://rdrr.io/r/base/readLines.html) on uncompressed
  humann mapping files.

## Value

a two-column `data.frame` that can be converted into an adjacency matrix
used as input for
[`weaveWeb()`](https://thomazbastiaanssen.github.io/anansi/reference/weaveWeb-generic.md).

## See also

[`weaveWeb()`](https://thomazbastiaanssen.github.io/anansi/reference/weaveWeb-generic.md)

## Examples

``` r
# some dummy input, as a character vector of IDs separated by tabs.
x <- c("x_1\ty_1\ty_2\ty_4", "x_2\ty_1\ty_3", "x_3\ty_2\ty_y")
linkBiobakeryMap(x)
#>   id.x id.y
#> 1  x_1  y_1
#> 2  x_1  y_2
#> 3  x_1  y_4
#> 4  x_2  y_1
#> 5  x_2  y_3
#> 6  x_3  y_2
#> 7  x_3  y_y
```
