# Methods for AnansiWeb S7 container class

Methods for AnansiWeb S7 container class

## Arguments

- x:

  input, AnansiWeb object

## Value

The desired information from an AnansiWeb object

## See also

- [`weaveWeb()`](https://thomazbastiaanssen.github.io/anansi/reference/weaveWeb-generic.md):
  for general use.

- [AnansiWeb-pairwise](https://thomazbastiaanssen.github.io/anansi/reference/AnansiWeb-pairwise.md):
  for methods for pairwise operations

## Examples

``` r
# Setup
web <- randomWeb(n_samp = 36)

# Accessors
dimnames(web)
#> $y
#>  [1] "y_1"  "y_2"  "y_3"  "y_4"  "y_5"  "y_6"  "y_7"  "y_8"  "y_9"  "y_10"
#> [11] "y_11" "y_12"
#> 
#> $x
#> [1] "x_1" "x_2" "x_3" "x_4" "x_5" "x_6" "x_7" "x_8"
#> 
dim(web)
#> [1] 12  8
names(web)
#> [1] "y" "x"

# Getters and setters:

tableX(web)[1:5, 1:5]
#>                       x
#> sample_id                      x_1        x_2       x_3         x_4         x_5
#>   anansi_ID_sample_1_1  0.88986536  0.1286881 0.6321729 -2.09497131 -1.74968202
#>   anansi_ID_sample_2_1 -1.48281332 -1.5330655 1.0849280  0.03407924 -0.03866763
#>   anansi_ID_sample_3_1  0.44575035  0.2023607 1.3564594  0.85272870 -0.79715324
#>   anansi_ID_sample_4_1  1.36977586 -0.7175387 0.3624240  0.74322814 -0.91592054
#>   anansi_ID_sample_5_1 -0.02011003  0.3616948 2.1693445  0.55715361  0.07326746
tableY(web)[1:5, 1:5]
#>                       y
#> sample_id                     y_1        y_2        y_3        y_4        y_5
#>   anansi_ID_sample_1_1  1.4154123  2.1264445 -0.3438315 -1.3416861  0.5662016
#>   anansi_ID_sample_2_1 -0.3837330 -1.4761969  1.0628765  1.1589792  1.1522120
#>   anansi_ID_sample_3_1 -0.1740864  0.4078885  0.8130582 -0.2032090 -0.7561974
#>   anansi_ID_sample_4_1 -0.2217445  1.3939778  1.8034834 -0.3780286 -0.4892583
#>   anansi_ID_sample_5_1 -1.0095287  0.3602783 -0.1050687  1.7361110 -1.1660523
dictionary(web)
#> 12 x 8 sparse Matrix of class "ngCMatrix"
#>       x
#> y      x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8
#>   y_1    |   .   |   |   |   |   .   .
#>   y_2    |   .   .   |   .   .   .   .
#>   y_3    |   .   .   .   .   .   |   .
#>   y_4    |   |   .   .   .   |   |   .
#>   y_5    |   |   .   .   |   |   .   |
#>   y_6    .   |   .   .   .   .   .   |
#>   y_7    .   |   .   .   .   |   |   |
#>   y_8    .   |   .   .   |   |   |   |
#>   y_9    |   |   |   |   |   .   |   |
#>   y_10   |   |   |   .   |   |   |   .
#>   y_11   .   |   .   .   .   |   .   .
#>   y_12   |   .   |   .   |   .   .   |
head(metadata(web))
#>                      sample_id repeated group_ab subtype     score_a    score_b
#> anansi_ID_sample_1_1  sample_1    rep_1        a       x -0.08013185 -0.6086578
#> anansi_ID_sample_2_1  sample_2    rep_1        a       y -0.03228341 -0.7311029
#> anansi_ID_sample_3_1  sample_3    rep_1        b       y -0.71898093  2.7151442
#> anansi_ID_sample_4_1  sample_4    rep_1        a       x -1.11656132 -1.3393870
#> anansi_ID_sample_5_1  sample_5    rep_1        a       z -0.78026990 -0.6460152
#> anansi_ID_sample_6_1  sample_6    rep_1        a       z -1.77695853 -0.9324546
#>                         score_c
#> anansi_ID_sample_1_1 -0.7432999
#> anansi_ID_sample_2_1 -0.3042238
#> anansi_ID_sample_3_1  0.3376581
#> anansi_ID_sample_4_1 -0.6075021
#> anansi_ID_sample_5_1 -0.2955603
#> anansi_ID_sample_6_1 -0.1345371

# Assign some random metadata
metadata(web) <- data.frame(
    id = row.names(tableY(web)),
    a = rnorm(36),
    b = sample(c("a", "b"), 36, TRUE),
    row.names = "id"
)

# Coerce to list
weblist <- as.list(web)

# Coerce to Data.frame
webdf <- as.data.frame(web)

# Coerce to MultiAssayExperiment
mae <- asMAE(web)

# Coerce to TreeSummarizedExperiment
tse <- asTSE(web)
```
