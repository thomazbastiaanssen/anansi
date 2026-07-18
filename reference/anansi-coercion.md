# Coercion functions for anansi

Coercion functions for anansi

## Usage

``` r
# S3 method for class '`anansi::AnansiWeb`'
as.list(x, ...)

# S3 method for class '`anansi::MultiFactor`'
as.list(x, ..., use.names = TRUE)

# S3 method for class '`anansi::AnansiWeb`'
as.data.frame(x, row.names, optional, ...)

asMAE(x)

asTSE(x)
```

## Arguments

- x:

  input object

- ...:

  additional arguments (currently not used).

- use.names:

  `Logical scalar`, whether output list should contain character
  (Default) or integer data frame. If `FALSE`, returns `unfactor(x)`.

- row.names, optional:

  Not used. See ?base::data.frame

## Value

An object of the desired class.

## See also

[`unfactor()`](https://thomazbastiaanssen.github.io/anansi/reference/MultiFactor-methods.md)

## Examples

``` r
# AnansiWeb
x <- randomWeb(
    n_samples = 5,
    n_features_x = 4,
    n_features_y = 6
)

as.list(x)
#> $y
#>                       y
#> sample_id                     y_1         y_2        y_3        y_4         y_5
#>   anansi_ID_sample_1_1 -0.4465888  0.28961148 -0.4051922 -0.5404371  2.58014889
#>   anansi_ID_sample_2_1  0.3371416 -0.13567761 -0.6837008  1.0204631  1.80346755
#>   anansi_ID_sample_3_1  1.0083097 -0.15330647  0.7336743 -1.0043337 -1.23211530
#>   anansi_ID_sample_4_1  2.4333666 -1.01651397 -1.6875778 -0.1178059 -0.06505926
#>   anansi_ID_sample_5_1  0.6852665 -0.06133107  0.3343705  0.4952910 -1.13000150
#>                       y
#> sample_id                      y_6
#>   anansi_ID_sample_1_1  0.82501143
#>   anansi_ID_sample_2_1  0.59267509
#>   anansi_ID_sample_3_1 -0.66034799
#>   anansi_ID_sample_4_1  0.05987437
#>   anansi_ID_sample_5_1 -1.05545080
#> 
#> $x
#>                       x
#> sample_id                     x_1        x_2        x_3        x_4
#>   anansi_ID_sample_1_1  0.9882897 -1.0876967 -0.1905976  0.6162500
#>   anansi_ID_sample_2_1 -1.3504865 -1.3436343 -0.5037913 -0.4533584
#>   anansi_ID_sample_3_1 -1.1244198 -0.3595517 -0.4692892  0.5229102
#>   anansi_ID_sample_4_1  0.2299108 -1.2574566  1.0016837  0.3534007
#>   anansi_ID_sample_5_1 -1.0237417 -1.1576458  1.2734635  0.2251215
#> 
#> $dictionary
#> 6 x 4 sparse Matrix of class "ngCMatrix"
#>      x
#> y     x_1 x_2 x_3 x_4
#>   y_1   |   .   .   |
#>   y_2   |   .   |   .
#>   y_3   .   .   |   |
#>   y_4   |   |   |   |
#>   y_5   |   .   .   .
#>   y_6   .   .   .   |
#> 
#> $metadata
#>                      sample_id repeated group_ab subtype    score_a    score_b
#> anansi_ID_sample_1_1  sample_1    rep_1        b       x -0.9799220 -0.7084848
#> anansi_ID_sample_2_1  sample_2    rep_1        a       x -1.4009123  1.4750925
#> anansi_ID_sample_3_1  sample_3    rep_1        a       y  1.4450149  0.8450042
#> anansi_ID_sample_4_1  sample_4    rep_1        b       x -0.4234816  1.2939944
#> anansi_ID_sample_5_1  sample_5    rep_1        a       y -0.7337496  0.2981611
#>                          score_c
#> anansi_ID_sample_1_1 -0.40568249
#> anansi_ID_sample_2_1 -0.13880758
#> anansi_ID_sample_3_1 -0.22259273
#> anansi_ID_sample_4_1  1.74715026
#> anansi_ID_sample_5_1 -0.08170467
#> 
as.data.frame(x)
#>                             y_1         y_2        y_3        y_4         y_5
#> anansi_ID_sample_1_1 -0.4465888  0.28961148 -0.4051922 -0.5404371  2.58014889
#> anansi_ID_sample_2_1  0.3371416 -0.13567761 -0.6837008  1.0204631  1.80346755
#> anansi_ID_sample_3_1  1.0083097 -0.15330647  0.7336743 -1.0043337 -1.23211530
#> anansi_ID_sample_4_1  2.4333666 -1.01651397 -1.6875778 -0.1178059 -0.06505926
#> anansi_ID_sample_5_1  0.6852665 -0.06133107  0.3343705  0.4952910 -1.13000150
#>                              y_6        x_1        x_2        x_3        x_4
#> anansi_ID_sample_1_1  0.82501143  0.9882897 -1.0876967 -0.1905976  0.6162500
#> anansi_ID_sample_2_1  0.59267509 -1.3504865 -1.3436343 -0.5037913 -0.4533584
#> anansi_ID_sample_3_1 -0.66034799 -1.1244198 -0.3595517 -0.4692892  0.5229102
#> anansi_ID_sample_4_1  0.05987437  0.2299108 -1.2574566  1.0016837  0.3534007
#> anansi_ID_sample_5_1 -1.05545080 -1.0237417 -1.1576458  1.2734635  0.2251215
#>                      sample_id repeated group_ab subtype    score_a    score_b
#> anansi_ID_sample_1_1  sample_1    rep_1        b       x -0.9799220 -0.7084848
#> anansi_ID_sample_2_1  sample_2    rep_1        a       x -1.4009123  1.4750925
#> anansi_ID_sample_3_1  sample_3    rep_1        a       y  1.4450149  0.8450042
#> anansi_ID_sample_4_1  sample_4    rep_1        b       x -0.4234816  1.2939944
#> anansi_ID_sample_5_1  sample_5    rep_1        a       y -0.7337496  0.2981611
#>                          score_c
#> anansi_ID_sample_1_1 -0.40568249
#> anansi_ID_sample_2_1 -0.13880758
#> anansi_ID_sample_3_1 -0.22259273
#> anansi_ID_sample_4_1  1.74715026
#> anansi_ID_sample_5_1 -0.08170467

# AnansiWeb to MultiAssayExperiment
asMAE(x)
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] y: SummarizedExperiment with 6 rows and 5 columns
#>  [2] x: SummarizedExperiment with 4 rows and 5 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
asTSE(x)
#> class: TreeSummarizedExperiment 
#> dim: 6 5 
#> metadata(1): dictionary
#> assays(1): y
#> rownames(6): y_1 y_2 ... y_5 y_6
#> rowData names(0):
#> colnames(5): anansi_ID_sample_1_1 anansi_ID_sample_2_1
#>   anansi_ID_sample_3_1 anansi_ID_sample_4_1 anansi_ID_sample_5_1
#> colData names(7): sample_id repeated ... score_b score_c
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(1): x
#> rowLinks: NULL
#> rowTree: NULL
#> colLinks: NULL
#> colTree: NULL

# MultiFactor
x <- randomMultiFactor(n_types = 3, n_features = 3)
as.list(x, use.names = TRUE)
#> $a2b
#>       a     b
#> 1 a_001 b_002
#> 2 a_002 b_002
#> 3 a_003 b_002
#> 4 a_001 b_003
#> 5 a_003 b_003
#> 
#> $b2c
#>       b     c
#> 1 b_002 c_001
#> 2 b_003 c_002
#> 3 b_002 c_003
#> 4 b_003 c_003
#> 
```
