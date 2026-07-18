# AnansiWeb S7 container class

`AnansiWeb` is an S7 class containing two feature tables as well as a
dictionary to link them. `AnansiWeb` is the main container that will
hold your input data throughout the `anansi` pipeline.

Typical use of the `anansi` package will involve generating an
`AnansiWeb` object using the
[`weaveWeb()`](https://thomazbastiaanssen.github.io/anansi/reference/weaveWeb-generic.md)
function.

The function `AnansiWeb()` constructs an `AnansiWeb` object from two
feature tables and an adjacency matrix.

## Usage

``` r
## Constructor for `AnansiWeb` objects
AnansiWeb(tableX, tableY, dictionary, metadata = data.frame())
```

## Arguments

- tableY, tableX:

  A table containing features of interest. Rows should be samples and
  columns should be features. Y and X refer to the position of the
  features in a formula: Y ~ X.

- dictionary:

  A binary adjacency matrix of class `Matrix`, or coercible to `Matrix`

- metadata:

  `list` of metadata. Optional.

## Value

an `AnansiWeb` object, with sparse binary biadjacency matrix with
features from `y` as rows and features from `x` as columns in
`dictionary` slot.

## Slots

- `tableY,tableX`:

  Two `matrix` objects of measurements, data. Rows are samples and
  columns are features. Access with `@tableY` and `@tableX`.

- `dictionary`:

  `Matrix`, binary adjacency matrix. Optionally sparse. Typically
  generated using the
  [`weaveWeb()`](https://thomazbastiaanssen.github.io/anansi/reference/weaveWeb-generic.md)
  function. Access with `@dictionary`.

- `metadata`:

  Optional `data.frame` of sample metadata. Access with `@metadata`.

## See also

- [AnansiWeb-methods](https://thomazbastiaanssen.github.io/anansi/reference/AnansiWeb-methods.md)

- [randomAnansi](https://thomazbastiaanssen.github.io/anansi/reference/randomAnansi.md):
  For generation of random AnansiWeb objects.

- [`kegg_link()`](https://thomazbastiaanssen.github.io/anansi/reference/kegg_link.md):
  For examples of input for link argument.

## Examples

``` r

# Use AnansiWeb() to construct an AnansiWeb object from components:
tX <- `dimnames<-`(replicate(5, (rnorm(36))),
    value = list(
        as.character(seq_len(36)),
        letters[1:5]
    )
)
tY <- `dimnames<-`(replicate(3, (rnorm(36))),
    value = list(
        as.character(seq_len(36)),
        LETTERS[1:3]
    )
)

d <- matrix(TRUE,
    nrow = NCOL(tY), ncol = NCOL(tX),

    # Note: Dictionary should have named dimensions
    dimnames = list(
        y_names = colnames(tY),
        x_names = colnames(tX)
    )
)
web <- AnansiWeb(tableX = tX, tableY = tY, dictionary = d)
#> Warning: Argument `metadata` not provided; Please validate sample ID order.
```
