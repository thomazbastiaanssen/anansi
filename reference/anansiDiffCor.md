# Run differential correlation analysis for all interacting metabolites and functions.

Can either take continuous or categorical data. Typically, the main
[`anansi()`](https://thomazbastiaanssen.github.io/anansi/reference/anansi.md)
function will run this for you.

## Usage

``` r
anansiDiffCor(web, sat_model, errorterm, int.terms, metadata, verbose)
```

## Arguments

- web:

  An `AnansiWeb` object, containing two tables with omics data and a
  dictionary that links them. See `weaveWebFromTables()` for how to
  weave a web.

- sat_model:

  A `formula` object, containing the full model

- errorterm:

  A `character vector`, containing the metadata col.name denoting
  repeated measures.

- int.terms:

  A `character vector`, containing the metadata col.names denoting
  covariates interacting with X to be tested for differential
  associations.

- metadata:

  A vector or data.frame of categorical or continuous value necessary
  for differential correlations. Often a state or treatment score. If no
  argument provided, anansi will let you know and still to regular
  correlations according to your dictionary.

- verbose:

  A boolean. Toggles whether to print diagnostic information while
  running. Useful for debugging errors on large datasets.

## Value

a list of `AnansiTale` result objects, one for the total model, one for
emergent correlations and one for disjointed correlations.
