# Generate a random AnansiWeb or MultiFactor

Randomly generate a valid `AnansiWeb` or `MultiFactor` object.

## Usage

``` r
randomWeb(
  n_samples = 10,
  n_reps = 1L,
  n_features_x = 8,
  n_features_y = 12,
  sparseness = 0.5,
  tableY = NULL,
  tableX = NULL,
  dictionary = NULL
)

randomMultiFactor(n_types = 6, n_features = 100, sparseness = 0.5)

krebsDemoWeb(n_samples = 100, n_reps = 4L)
```

## Arguments

- n_samples, n_reps:

  `Numeric scalar` Number of samples and repeated measures of those
  samples to be generated. Ignored if `tableY` and `tableX` are
  provided. (defaults: 10 samples without repeats)

- n_features_y, n_features_x:

  `Numeric scalar` Number of features to be generated. Ignored if
  `tableY` and `tableX` are provided.

- sparseness:

  `Numeric scalar`, proportion: How rare are connections

- tableY, tableX:

  A table containing features of interest. Rows should be samples and
  columns should be features. Y and X refer to the position of the
  features in a formula: Y ~ X.

- dictionary:

  A binary adjacency matrix of class `Matrix`, or coercible to `Matrix`.

- n_types:

  `Numeric scalar`, number of types of features to generate

- n_features:

  `Numeric scalar`, number of features per type

## Value

a randomly generated object of the specified class.

## See also

[`AnansiWeb()`](https://thomazbastiaanssen.github.io/anansi/reference/AnansiWeb-class.md),
[`MultiFactor()`](https://thomazbastiaanssen.github.io/anansi/reference/MultiFactor-class.md)

## Examples

``` r
# Make a random AnansiWeb object
randomWeb()
#> anansi::AnansiWeb S7_object with 10 observations:
#>     tableY: y (12 features)
#>     tableX: x (8 features)
#> Use @ to access: tableX, tableY, dictionary, metadata.
krebsDemoWeb()
#> anansi::AnansiWeb S7_object with 400 observations:
#>     tableY: Enzyme (8 features)
#>     tableX: Metabolite (9 features)
#> Use @ to access: tableX, tableY, dictionary, metadata.
randomMultiFactor()
#> An anansi::MultiFactor S7_object,
#>     6 feature types across 5 edge lists.
#> 
#>       a   b   c   d   e   f
#> a2b 100 100   .   .   .   .
#> b2c   . 100 100   .   .   .
#> c2d   .   . 100 100   .   .
#> d2e   .   .   . 100 100   .
#> e2f   .   .   .   . 100 100
#> 
#> Values represent unique feature names in that edge list.
#> 
#> Levels:
#> 
#> a : 100 Levels: a_001 a_002 ... a_100 
#> b : 100 Levels: b_001 b_002 ... b_100 
#> c : 100 Levels: c_001 c_002 ... c_100 
#> d : 100 Levels: d_001 d_002 ... d_100 
#> e : 100 Levels: e_001 e_002 ... e_100 
#> f : 100 Levels: f_001 f_002 ... f_100 
```
