# Get a list of all pairs of features

Get a list of all pairs of features

## Usage

``` r
getFeaturePairs(x, ...)
```

## Arguments

- x:

  input object

- ...:

  additional arguments for specific methods

## Value

an list of two-column data.frames that represent all feature pairs.

## Examples

``` r
x <- randomWeb(10)
head(getFeaturePairs(x), 3)
#> [[1]]
#>                              y_1        x_1
#> anansi_ID_sample_1_1   1.6787554  0.8532586
#> anansi_ID_sample_2_1  -0.2514572 -0.8401742
#> anansi_ID_sample_3_1  -0.8552489 -0.2972984
#> anansi_ID_sample_4_1   0.0293183  1.0433654
#> anansi_ID_sample_5_1  -0.3520895 -2.6550646
#> anansi_ID_sample_6_1  -1.3680452 -0.1460115
#> anansi_ID_sample_7_1   0.4846134  0.1121593
#> anansi_ID_sample_8_1  -0.7003770  0.3439378
#> anansi_ID_sample_9_1  -0.1648789  1.2877091
#> anansi_ID_sample_10_1  0.4580152  0.8977251
#> 
#> [[2]]
#>                              y_6        x_1
#> anansi_ID_sample_1_1   1.3764592  0.8532586
#> anansi_ID_sample_2_1   0.7078778 -0.8401742
#> anansi_ID_sample_3_1  -2.0243059 -0.2972984
#> anansi_ID_sample_4_1  -0.8033473  1.0433654
#> anansi_ID_sample_5_1  -0.1740919 -2.6550646
#> anansi_ID_sample_6_1   0.5078751 -0.1460115
#> anansi_ID_sample_7_1   0.5992470  0.1121593
#> anansi_ID_sample_8_1   1.4517137  0.3439378
#> anansi_ID_sample_9_1   1.3859527  1.2877091
#> anansi_ID_sample_10_1  0.9360747  0.8977251
#> 
#> [[3]]
#>                              y_10        x_1
#> anansi_ID_sample_1_1  -0.02622092  0.8532586
#> anansi_ID_sample_2_1  -1.13863372 -0.8401742
#> anansi_ID_sample_3_1   0.68138608 -0.2972984
#> anansi_ID_sample_4_1  -0.84193128  1.0433654
#> anansi_ID_sample_5_1  -0.50130139 -2.6550646
#> anansi_ID_sample_6_1  -1.25781113 -0.1460115
#> anansi_ID_sample_7_1   2.08118340  0.1121593
#> anansi_ID_sample_8_1  -0.99627177  0.3439378
#> anansi_ID_sample_9_1  -1.34681825  1.2877091
#> anansi_ID_sample_10_1 -0.92133934  0.8977251
#> 
```
