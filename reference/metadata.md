# Get metadata.

Get metadata.

## Arguments

- x:

  input object

- ...:

  additional arguments

## Value

`metadata` slot of x

## Details

Compatible with S4Vectors generic.

## Examples

``` r
x <- randomWeb(10)
metadata(x)
#>                       sample_id repeated group_ab subtype     score_a
#> anansi_ID_sample_1_1   sample_1    rep_1        a       x  0.89176234
#> anansi_ID_sample_2_1   sample_2    rep_1        a       x -0.13629434
#> anansi_ID_sample_3_1   sample_3    rep_1        a       z  0.28119119
#> anansi_ID_sample_4_1   sample_4    rep_1        a       z -0.11081246
#> anansi_ID_sample_5_1   sample_5    rep_1        a       z -0.03624743
#> anansi_ID_sample_6_1   sample_6    rep_1        b       y  0.02497722
#> anansi_ID_sample_7_1   sample_7    rep_1        a       x  0.13105356
#> anansi_ID_sample_8_1   sample_8    rep_1        b       y  2.59552392
#> anansi_ID_sample_9_1   sample_9    rep_1        a       z  0.41928399
#> anansi_ID_sample_10_1 sample_10    rep_1        a       z -0.28169842
#>                           score_b     score_c
#> anansi_ID_sample_1_1   0.54486311  0.10347211
#> anansi_ID_sample_2_1   0.86823069 -0.65763660
#> anansi_ID_sample_3_1  -0.20667527  0.03387470
#> anansi_ID_sample_4_1   0.07494732 -0.64975656
#> anansi_ID_sample_5_1  -2.25090133  0.91123965
#> anansi_ID_sample_6_1   1.50096116 -0.04727292
#> anansi_ID_sample_7_1  -0.48276894 -1.17851501
#> anansi_ID_sample_8_1   0.60521064  2.24266097
#> anansi_ID_sample_9_1   0.80068548  1.46748561
#> anansi_ID_sample_10_1  0.31836219  0.55206173
```
