# Weave an AnansiWeb object

Weave an AnansiWeb object

## Usage

``` r
weaveWeb(x, ...)

weaveKEGG(x, ...)
```

## Arguments

- x:

  input object

- ...:

  additional arguments

## Value

an `AnansiWeb` object, with sparse binary biadjacency matrix with
features from `y` as rows and features from `x` as columns in
`dictionary` slot.

## See also

[`weaveWeb-methods()`](https://thomazbastiaanssen.github.io/anansi/reference/weaveWeb-methods.md)

## Examples

``` r
# Setup demo tables
ec2ko <- kegg_link()[["ec2ko"]]
ec2cpd <- kegg_link()[["ec2cpd"]]

# Basic usage
weaveWeb(cpd ~ ko, link = kegg_link())
#> Warning: Argument `metadata` not provided; Please validate sample ID order.
#> anansi::AnansiWeb S7_object with 1 observations:
#>     tableY: cpd (6754 features)
#>     tableX: ko (7270 features)
#> Use @ to access: tableX, tableY, dictionary, metadata.
weaveWeb(x = "ko", y = "ec", link = ec2ko)
#> Warning: Argument `metadata` not provided; Please validate sample ID order.
#> anansi::AnansiWeb S7_object with 1 observations:
#>     tableY: ec (5581 features)
#>     tableX: ko (8817 features)
#> Use @ to access: tableX, tableY, dictionary, metadata.
weaveWeb(ec ~ cpd, link = ec2cpd)
#> Warning: Argument `metadata` not provided; Please validate sample ID order.
#> anansi::AnansiWeb S7_object with 1 observations:
#>     tableY: ec (5934 features)
#>     tableX: cpd (7802 features)
#> Use @ to access: tableX, tableY, dictionary, metadata.

# A wrapper is available for kegg ko, ec and cpd data
generic <- weaveWeb(cpd ~ ko, link = kegg_link())
#> Warning: Argument `metadata` not provided; Please validate sample ID order.
kegg_wrapper <- weaveKEGG(cpd ~ ko)
#> Warning: Argument `metadata` not provided; Please validate sample ID order.

identical(generic, kegg_wrapper)
#> [1] TRUE

# The following are equivalent to transposition:
a <- weaveWeb(ko ~ cpd, link = kegg_link()) |> dictionary()
#> Warning: Argument `metadata` not provided; Please validate sample ID order.
b <- weaveWeb(cpd ~ ko, link = kegg_link()) |> dictionary()
#> Warning: Argument `metadata` not provided; Please validate sample ID order.

identical(a, Matrix::t(b))
#> [1] TRUE
```
