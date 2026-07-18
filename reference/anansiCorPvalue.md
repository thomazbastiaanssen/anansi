# Compute r-statistics for each featureY-featureX pair in the dictionary. Typically, the main `anansi()` function will run this for you.

Compute r-statistics for each featureY-featureX pair in the dictionary.
Typically, the main
[`anansi()`](https://thomazbastiaanssen.github.io/anansi/reference/anansi.md)
function will run this for you.

## Usage

``` r
anansiCorPvalue(web, group.bool, verbose)
```

## Arguments

- web:

  An `AnansiWeb` object, containing two tables with omics data and a
  dictionary that links them. See `weaveWebFromTables()` for how to
  weave a web.

- group.bool:

  A categorical or continuous value necessary for differential
  correlations. Typically a state or treatment score. If no argument
  provided, anansi will let you know and still to regular correlations
  according to your dictionary.

- verbose:

  A boolean. Toggles whether to print diagnostic information while
  running. Useful for debugging errors on large datasets.

## Value

An `AnansiTale` result object.

## See also

[`anansi()`](https://thomazbastiaanssen.github.io/anansi/reference/anansi.md)\
[`anansiCorTestByGroup()`](https://thomazbastiaanssen.github.io/anansi/reference/anansiCorTestByGroup.md)
