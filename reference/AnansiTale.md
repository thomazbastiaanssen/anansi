# AnansiTale S7 container class. Not intended for general use.

`AnansiTale` is the main container that will hold your stats output data
coming out of the `anansi` pipeline.

## Usage

``` r
AnansiTale(
  subject = character(0),
  type = character(0),
  df = integer(0),
  estimates = integer(0),
  f.values = integer(0),
  t.values = integer(0),
  p.values = integer(0)
)
```

## Arguments

- subject:

  A character that describes the data that was queried.

- type:

  A character that describes type of parameter contained in the
  `estimates` slot. For example r.values for correlations or r.squared
  for models.

- df:

  a vector of length 2, containing df1 and df2 corresponding to the
  F-ratio considered.

- estimates:

  A matrix containing the estimates for the parameters named in the
  `type` slot.

- f.values:

  A matrix containing the f-values, for least-squares.

- t.values:

  A matrix containing the t-values, for correlations.

- p.values:

  A matrix containing the p.values for the parameters named in the
  `type` slot.

## Value

an `AnansiTale`

## Slots

- `subject`:

  A character that describes the data that was queried.

- `type`:

  A character that describes type of parameter contained in the
  `estimates` slot. For example r.values for correlations or r.squared
  for models.

- `df`:

  a vector of length 2, containing df1 and df2 corresponding to the
  F-ratio considered.

- `estimates`:

  A matrix containing the estimates for the parameters named in the
  `type` slot.

- `f.values`:

  A matrix containing the f-values, for least-squares.

- `t.values`:

  A matrix containing the t-values, for correlations.

- `p.values`:

  A matrix containing the p.values for the parameters named in the
  `type` slot.

## Examples

``` r
AnansiTale
#> <anansi::AnansiTale> class
#> @ parent     : <S7_object>
#> @ constructor: function(subject, type, df, estimates, f.values, t.values, p.values) {...}
#> @ validator  : <NULL>
#> @ properties :
#>  $ subject  : <character>          
#>  $ type     : <character>          
#>  $ df       : <integer> or <double>
#>  $ estimates: <integer> or <double>
#>  $ f.values : <integer> or <double>
#>  $ t.values : <integer> or <double>
#>  $ p.values : <integer> or <double>
```
