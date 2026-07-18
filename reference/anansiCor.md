# Compute r-statistics for each featureY-featureX pair in the dictionary. Typically, the main `anansi()` function will run this for you.

Compute r-statistics for each featureY-featureX pair in the dictionary.
Typically, the main
[`anansi()`](https://thomazbastiaanssen.github.io/anansi/reference/anansi.md)
function will run this for you.

## Usage

``` r
anansiCor(web, group.bool)
```

## Arguments

- web:

  An `AnansiWeb` object, containing two tables with omics data and a
  dictionary that links them. See `weaveWebFromTables()` for how to
  weave a web.

- group.bool:

  A boolean vector used to select which samples should be included in
  the correlations.

## Value

A matrix of r-statistics.

## See also

[`anansi()`](https://thomazbastiaanssen.github.io/anansi/reference/anansi.md)\
[`anansiCorTestByGroup()`](https://thomazbastiaanssen.github.io/anansi/reference/anansiCorTestByGroup.md)
