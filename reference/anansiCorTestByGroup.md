# Run correlations for all interacting metabolites and functions.

If the `groups` argument is suitable, will also run correlation analysis
per group. Typically, the main
[`anansi()`](https://thomazbastiaanssen.github.io/anansi/reference/anansi.md)
function will run this for you.

## Usage

``` r
anansiCorTestByGroup(web, group.vec, verbose = TRUE)
```

## Arguments

- web:

  An `AnansiWeb` object, containing two tables with omics data and a
  dictionary that links them. See `weaveWebFromTables()` for how to
  weave a web.

- group.vec:

  A character vector denoting group membership. Typically a state or
  treatment score.

- verbose:

  A boolean. Toggles whether to print diagnostic information while
  running. Useful for debugging errors on large datasets.

## Value

a list of `AnansiTale` result objects, one for the total dataset and per
group if applicable.

## See also

[`anansi()`](https://thomazbastiaanssen.github.io/anansi/reference/anansi.md)
