# Set metadata.

Set metadata.

## Arguments

- x:

  input object

- ...:

  additional arguments

- value:

  replacement value, coerced to data.frame

## Value

`x` with modified `metadata` slot.

## Details

Compatible with S4Vectors generic.

## Examples

``` r
x <- randomWeb(10)
metadata(x) <- cbind(
    metadata(x),
    new_groups = c("A", "B")
)
```
