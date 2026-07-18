# Methods for pairwise operations on AnansiWeb objects

Methods to run pairwise analysis on multi-modal data.

## Arguments

- X, x:

  input, AnansiWeb object

- ...:

  additional arguments

- FUN:

  a function with at least two arguments. The variables `x` and `y`, in
  order, refer to the corresponding values of feature pairs in `tableX`
  and `tableY`.

- MoreArgs, SIMPLIFY, USE.NAMES:

  see ?base::mapply

- which:

  `integer matrix`, indicating pair positions in `x@tableY` and
  `x@tableX`, respectively. If `NULL` (default):
  `Matrix::which(x@dictionary, TRUE)`.

- with.metadata:

  `Logical scalar` whether to append metadata to output

## Value

pairs: A two-column array index, corresponding to i,j coordinates in
matrix notation.\
getFeaturePairs: A list of data.frames with the paired data\

## Examples

``` r
web <- randomWeb()
# Extract data.frames in pairs (only show first)
getFeaturePairs(web)[1L]
#> [[1]]
#>                              y_1        x_1
#> anansi_ID_sample_1_1   0.8250203  1.2122810
#> anansi_ID_sample_2_1  -0.7501828 -1.2511968
#> anansi_ID_sample_3_1   0.5362105  0.4002660
#> anansi_ID_sample_4_1  -1.5727175  0.9347076
#> anansi_ID_sample_5_1  -0.9864274 -0.4726698
#> anansi_ID_sample_6_1   1.9835123  1.5830132
#> anansi_ID_sample_7_1  -1.8513784  0.5977542
#> anansi_ID_sample_8_1  -0.9097882  0.4660846
#> anansi_ID_sample_9_1  -1.9514400 -0.7443675
#> anansi_ID_sample_10_1 -0.8002765  0.9590553
#> 

pairs(x = web)
#>       [,1] [,2]
#>  [1,]    1    1
#>  [2,]    3    1
#>  [3,]    4    1
#>  [4,]    9    1
#>  [5,]   11    1
#>  [6,]   12    1
#>  [7,]    1    2
#>  [8,]    2    2
#>  [9,]    4    2
#> [10,]    5    2
#> [11,]    8    2
#> [12,]    9    2
#> [13,]   11    2
#> [14,]    2    3
#> [15,]    4    3
#> [16,]    6    3
#> [17,]    7    3
#> [18,]    8    3
#> [19,]   10    3
#> [20,]   12    3
#> [21,]    1    4
#> [22,]    2    4
#> [23,]    3    4
#> [24,]    6    4
#> [25,]    9    4
#> [26,]   12    4
#> [27,]    6    5
#> [28,]    7    5
#> [29,]    9    5
#> [30,]   10    5
#> [31,]   12    5
#> [32,]    1    6
#> [33,]    2    6
#> [34,]    4    6
#> [35,]    5    6
#> [36,]    9    6
#> [37,]   11    6
#> [38,]   12    6
#> [39,]    1    7
#> [40,]    2    7
#> [41,]    6    7
#> [42,]    7    7
#> [43,]   11    7
#> [44,]   12    7
#> [45,]    1    8
#> [46,]    3    8
#> [47,]    4    8
#> [48,]    9    8

getFeaturePairs(x = web, which = NULL, with.metadata = FALSE)
#> [[1]]
#>                              y_1        x_1
#> anansi_ID_sample_1_1   0.8250203  1.2122810
#> anansi_ID_sample_2_1  -0.7501828 -1.2511968
#> anansi_ID_sample_3_1   0.5362105  0.4002660
#> anansi_ID_sample_4_1  -1.5727175  0.9347076
#> anansi_ID_sample_5_1  -0.9864274 -0.4726698
#> anansi_ID_sample_6_1   1.9835123  1.5830132
#> anansi_ID_sample_7_1  -1.8513784  0.5977542
#> anansi_ID_sample_8_1  -0.9097882  0.4660846
#> anansi_ID_sample_9_1  -1.9514400 -0.7443675
#> anansi_ID_sample_10_1 -0.8002765  0.9590553
#> 
#> [[2]]
#>                              y_3        x_1
#> anansi_ID_sample_1_1  -1.7368783  1.2122810
#> anansi_ID_sample_2_1   0.0404813 -1.2511968
#> anansi_ID_sample_3_1  -0.1241795  0.4002660
#> anansi_ID_sample_4_1  -0.6127534  0.9347076
#> anansi_ID_sample_5_1   0.1606625 -0.4726698
#> anansi_ID_sample_6_1  -0.6623596  1.5830132
#> anansi_ID_sample_7_1  -0.3348517  0.5977542
#> anansi_ID_sample_8_1   0.6230116  0.4660846
#> anansi_ID_sample_9_1   1.0281943 -0.7443675
#> anansi_ID_sample_10_1 -1.1345782  0.9590553
#> 
#> [[3]]
#>                              y_4        x_1
#> anansi_ID_sample_1_1   0.9169111  1.2122810
#> anansi_ID_sample_2_1   1.2138748 -1.2511968
#> anansi_ID_sample_3_1  -0.6859325  0.4002660
#> anansi_ID_sample_4_1  -1.0945682  0.9347076
#> anansi_ID_sample_5_1   0.3648741 -0.4726698
#> anansi_ID_sample_6_1  -1.0080171  1.5830132
#> anansi_ID_sample_7_1   0.5549964  0.5977542
#> anansi_ID_sample_8_1   0.4360659  0.4660846
#> anansi_ID_sample_9_1   0.1053733 -0.7443675
#> anansi_ID_sample_10_1 -0.2537379  0.9590553
#> 
#> [[4]]
#>                               y_9        x_1
#> anansi_ID_sample_1_1  -1.69018029  1.2122810
#> anansi_ID_sample_2_1  -0.94549766 -1.2511968
#> anansi_ID_sample_3_1  -0.34301196  0.4002660
#> anansi_ID_sample_4_1   0.35846241  0.9347076
#> anansi_ID_sample_5_1   0.04818886 -0.4726698
#> anansi_ID_sample_6_1   1.13017228  1.5830132
#> anansi_ID_sample_7_1  -0.52763231  0.5977542
#> anansi_ID_sample_8_1   0.47957533  0.4660846
#> anansi_ID_sample_9_1   1.41752930 -0.7443675
#> anansi_ID_sample_10_1 -0.08009077  0.9590553
#> 
#> [[5]]
#>                             y_11        x_1
#> anansi_ID_sample_1_1  -1.1431376  1.2122810
#> anansi_ID_sample_2_1   0.2334402 -1.2511968
#> anansi_ID_sample_3_1  -0.1892366  0.4002660
#> anansi_ID_sample_4_1  -1.6634113  0.9347076
#> anansi_ID_sample_5_1   1.9156908 -0.4726698
#> anansi_ID_sample_6_1  -0.8168944  1.5830132
#> anansi_ID_sample_7_1   0.3836505  0.5977542
#> anansi_ID_sample_8_1  -0.4586323  0.4660846
#> anansi_ID_sample_9_1  -0.7157040 -0.7443675
#> anansi_ID_sample_10_1  0.4295485  0.9590553
#> 
#> [[6]]
#>                              y_12        x_1
#> anansi_ID_sample_1_1  -1.16985274  1.2122810
#> anansi_ID_sample_2_1   1.32208831 -1.2511968
#> anansi_ID_sample_3_1   1.36596464  0.4002660
#> anansi_ID_sample_4_1  -0.24972371  0.9347076
#> anansi_ID_sample_5_1  -0.42472135 -0.4726698
#> anansi_ID_sample_6_1   1.09843934  1.5830132
#> anansi_ID_sample_7_1   0.67092304  0.5977542
#> anansi_ID_sample_8_1   0.03095466  0.4660846
#> anansi_ID_sample_9_1   2.36851196 -0.7443675
#> anansi_ID_sample_10_1  2.58578811  0.9590553
#> 
#> [[7]]
#>                              y_1        x_2
#> anansi_ID_sample_1_1   0.8250203  0.1526371
#> anansi_ID_sample_2_1  -0.7501828 -0.7693687
#> anansi_ID_sample_3_1   0.5362105  1.3930807
#> anansi_ID_sample_4_1  -1.5727175  0.9186688
#> anansi_ID_sample_5_1  -0.9864274  1.0160871
#> anansi_ID_sample_6_1   1.9835123  0.5346354
#> anansi_ID_sample_7_1  -1.8513784 -0.9876286
#> anansi_ID_sample_8_1  -0.9097882  1.1877923
#> anansi_ID_sample_9_1  -1.9514400 -0.5175422
#> anansi_ID_sample_10_1 -0.8002765 -0.2595603
#> 
#> [[8]]
#>                               y_2        x_2
#> anansi_ID_sample_1_1  -1.86936160  0.1526371
#> anansi_ID_sample_2_1  -0.75057148 -0.7693687
#> anansi_ID_sample_3_1  -0.59093839  1.3930807
#> anansi_ID_sample_4_1  -0.74209927  0.9186688
#> anansi_ID_sample_5_1   0.69351208  1.0160871
#> anansi_ID_sample_6_1  -0.05946397  0.5346354
#> anansi_ID_sample_7_1  -1.86389510 -0.9876286
#> anansi_ID_sample_8_1  -1.27450892  1.1877923
#> anansi_ID_sample_9_1  -1.78177081 -0.5175422
#> anansi_ID_sample_10_1 -0.50859339 -0.2595603
#> 
#> [[9]]
#>                              y_4        x_2
#> anansi_ID_sample_1_1   0.9169111  0.1526371
#> anansi_ID_sample_2_1   1.2138748 -0.7693687
#> anansi_ID_sample_3_1  -0.6859325  1.3930807
#> anansi_ID_sample_4_1  -1.0945682  0.9186688
#> anansi_ID_sample_5_1   0.3648741  1.0160871
#> anansi_ID_sample_6_1  -1.0080171  0.5346354
#> anansi_ID_sample_7_1   0.5549964 -0.9876286
#> anansi_ID_sample_8_1   0.4360659  1.1877923
#> anansi_ID_sample_9_1   0.1053733 -0.5175422
#> anansi_ID_sample_10_1 -0.2537379 -0.2595603
#> 
#> [[10]]
#>                              y_5        x_2
#> anansi_ID_sample_1_1  -0.9764912  0.1526371
#> anansi_ID_sample_2_1   0.3942178 -0.7693687
#> anansi_ID_sample_3_1   1.0149030  1.3930807
#> anansi_ID_sample_4_1   0.6314438  0.9186688
#> anansi_ID_sample_5_1  -0.4833983  1.0160871
#> anansi_ID_sample_6_1  -1.5263339  0.5346354
#> anansi_ID_sample_7_1   0.6168423 -0.9876286
#> anansi_ID_sample_8_1   1.2292868  1.1877923
#> anansi_ID_sample_9_1   0.1574530 -0.5175422
#> anansi_ID_sample_10_1  1.4063421 -0.2595603
#> 
#> [[11]]
#>                               y_8        x_2
#> anansi_ID_sample_1_1  -0.98534392  0.1526371
#> anansi_ID_sample_2_1  -1.49118975 -0.7693687
#> anansi_ID_sample_3_1   0.81959245  1.3930807
#> anansi_ID_sample_4_1   1.01340453  0.9186688
#> anansi_ID_sample_5_1   1.05360525  1.0160871
#> anansi_ID_sample_6_1  -0.07142864  0.5346354
#> anansi_ID_sample_7_1   0.95244461 -0.9876286
#> anansi_ID_sample_8_1  -0.38403752  1.1877923
#> anansi_ID_sample_9_1  -0.41775999 -0.5175422
#> anansi_ID_sample_10_1 -0.47026817 -0.2595603
#> 
#> [[12]]
#>                               y_9        x_2
#> anansi_ID_sample_1_1  -1.69018029  0.1526371
#> anansi_ID_sample_2_1  -0.94549766 -0.7693687
#> anansi_ID_sample_3_1  -0.34301196  1.3930807
#> anansi_ID_sample_4_1   0.35846241  0.9186688
#> anansi_ID_sample_5_1   0.04818886  1.0160871
#> anansi_ID_sample_6_1   1.13017228  0.5346354
#> anansi_ID_sample_7_1  -0.52763231 -0.9876286
#> anansi_ID_sample_8_1   0.47957533  1.1877923
#> anansi_ID_sample_9_1   1.41752930 -0.5175422
#> anansi_ID_sample_10_1 -0.08009077 -0.2595603
#> 
#> [[13]]
#>                             y_11        x_2
#> anansi_ID_sample_1_1  -1.1431376  0.1526371
#> anansi_ID_sample_2_1   0.2334402 -0.7693687
#> anansi_ID_sample_3_1  -0.1892366  1.3930807
#> anansi_ID_sample_4_1  -1.6634113  0.9186688
#> anansi_ID_sample_5_1   1.9156908  1.0160871
#> anansi_ID_sample_6_1  -0.8168944  0.5346354
#> anansi_ID_sample_7_1   0.3836505 -0.9876286
#> anansi_ID_sample_8_1  -0.4586323  1.1877923
#> anansi_ID_sample_9_1  -0.7157040 -0.5175422
#> anansi_ID_sample_10_1  0.4295485 -0.2595603
#> 
#> [[14]]
#>                               y_2         x_3
#> anansi_ID_sample_1_1  -1.86936160 -0.32806467
#> anansi_ID_sample_2_1  -0.75057148  0.07343239
#> anansi_ID_sample_3_1  -0.59093839 -0.24786302
#> anansi_ID_sample_4_1  -0.74209927 -1.37386226
#> anansi_ID_sample_5_1   0.69351208 -0.04044582
#> anansi_ID_sample_6_1  -0.05946397  0.42153824
#> anansi_ID_sample_7_1  -1.86389510  0.20159751
#> anansi_ID_sample_8_1  -1.27450892 -1.69719192
#> anansi_ID_sample_9_1  -1.78177081  0.64228768
#> anansi_ID_sample_10_1 -0.50859339 -0.99523961
#> 
#> [[15]]
#>                              y_4         x_3
#> anansi_ID_sample_1_1   0.9169111 -0.32806467
#> anansi_ID_sample_2_1   1.2138748  0.07343239
#> anansi_ID_sample_3_1  -0.6859325 -0.24786302
#> anansi_ID_sample_4_1  -1.0945682 -1.37386226
#> anansi_ID_sample_5_1   0.3648741 -0.04044582
#> anansi_ID_sample_6_1  -1.0080171  0.42153824
#> anansi_ID_sample_7_1   0.5549964  0.20159751
#> anansi_ID_sample_8_1   0.4360659 -1.69719192
#> anansi_ID_sample_9_1   0.1053733  0.64228768
#> anansi_ID_sample_10_1 -0.2537379 -0.99523961
#> 
#> [[16]]
#>                               y_6         x_3
#> anansi_ID_sample_1_1   0.39680092 -0.32806467
#> anansi_ID_sample_2_1   0.32642296  0.07343239
#> anansi_ID_sample_3_1  -1.45934331 -0.24786302
#> anansi_ID_sample_4_1  -0.79929553 -1.37386226
#> anansi_ID_sample_5_1   0.47001280 -0.04044582
#> anansi_ID_sample_6_1  -1.15018052  0.42153824
#> anansi_ID_sample_7_1  -0.27159672  0.20159751
#> anansi_ID_sample_8_1   0.45742422 -1.69719192
#> anansi_ID_sample_9_1  -0.01695794  0.64228768
#> anansi_ID_sample_10_1 -0.54159436 -0.99523961
#> 
#> [[17]]
#>                               y_7         x_3
#> anansi_ID_sample_1_1   0.87581342 -0.32806467
#> anansi_ID_sample_2_1   0.74083815  0.07343239
#> anansi_ID_sample_3_1  -0.16481239 -0.24786302
#> anansi_ID_sample_4_1  -0.75096428 -1.37386226
#> anansi_ID_sample_5_1  -1.25336629 -0.04044582
#> anansi_ID_sample_6_1  -1.03623626  0.42153824
#> anansi_ID_sample_7_1  -0.02943857  0.20159751
#> anansi_ID_sample_8_1  -0.27494378 -1.69719192
#> anansi_ID_sample_9_1   0.48545147  0.64228768
#> anansi_ID_sample_10_1  0.91698203 -0.99523961
#> 
#> [[18]]
#>                               y_8         x_3
#> anansi_ID_sample_1_1  -0.98534392 -0.32806467
#> anansi_ID_sample_2_1  -1.49118975  0.07343239
#> anansi_ID_sample_3_1   0.81959245 -0.24786302
#> anansi_ID_sample_4_1   1.01340453 -1.37386226
#> anansi_ID_sample_5_1   1.05360525 -0.04044582
#> anansi_ID_sample_6_1  -0.07142864  0.42153824
#> anansi_ID_sample_7_1   0.95244461  0.20159751
#> anansi_ID_sample_8_1  -0.38403752 -1.69719192
#> anansi_ID_sample_9_1  -0.41775999  0.64228768
#> anansi_ID_sample_10_1 -0.47026817 -0.99523961
#> 
#> [[19]]
#>                             y_10         x_3
#> anansi_ID_sample_1_1  -0.3941870 -0.32806467
#> anansi_ID_sample_2_1   0.5252576  0.07343239
#> anansi_ID_sample_3_1  -0.7959351 -0.24786302
#> anansi_ID_sample_4_1   0.2450964 -1.37386226
#> anansi_ID_sample_5_1  -0.7793198 -0.04044582
#> anansi_ID_sample_6_1   0.5949239  0.42153824
#> anansi_ID_sample_7_1   1.1088909  0.20159751
#> anansi_ID_sample_8_1  -0.9429690 -1.69719192
#> anansi_ID_sample_9_1   0.7003440  0.64228768
#> anansi_ID_sample_10_1 -0.4246748 -0.99523961
#> 
#> [[20]]
#>                              y_12         x_3
#> anansi_ID_sample_1_1  -1.16985274 -0.32806467
#> anansi_ID_sample_2_1   1.32208831  0.07343239
#> anansi_ID_sample_3_1   1.36596464 -0.24786302
#> anansi_ID_sample_4_1  -0.24972371 -1.37386226
#> anansi_ID_sample_5_1  -0.42472135 -0.04044582
#> anansi_ID_sample_6_1   1.09843934  0.42153824
#> anansi_ID_sample_7_1   0.67092304  0.20159751
#> anansi_ID_sample_8_1   0.03095466 -1.69719192
#> anansi_ID_sample_9_1   2.36851196  0.64228768
#> anansi_ID_sample_10_1  2.58578811 -0.99523961
#> 
#> [[21]]
#>                              y_1        x_4
#> anansi_ID_sample_1_1   0.8250203  0.9638139
#> anansi_ID_sample_2_1  -0.7501828 -1.6560372
#> anansi_ID_sample_3_1   0.5362105  1.0708611
#> anansi_ID_sample_4_1  -1.5727175 -0.1090264
#> anansi_ID_sample_5_1  -0.9864274  1.8991864
#> anansi_ID_sample_6_1   1.9835123 -1.1370307
#> anansi_ID_sample_7_1  -1.8513784 -0.2797198
#> anansi_ID_sample_8_1  -0.9097882 -0.8941291
#> anansi_ID_sample_9_1  -1.9514400  0.1367018
#> anansi_ID_sample_10_1 -0.8002765 -0.7491654
#> 
#> [[22]]
#>                               y_2        x_4
#> anansi_ID_sample_1_1  -1.86936160  0.9638139
#> anansi_ID_sample_2_1  -0.75057148 -1.6560372
#> anansi_ID_sample_3_1  -0.59093839  1.0708611
#> anansi_ID_sample_4_1  -0.74209927 -0.1090264
#> anansi_ID_sample_5_1   0.69351208  1.8991864
#> anansi_ID_sample_6_1  -0.05946397 -1.1370307
#> anansi_ID_sample_7_1  -1.86389510 -0.2797198
#> anansi_ID_sample_8_1  -1.27450892 -0.8941291
#> anansi_ID_sample_9_1  -1.78177081  0.1367018
#> anansi_ID_sample_10_1 -0.50859339 -0.7491654
#> 
#> [[23]]
#>                              y_3        x_4
#> anansi_ID_sample_1_1  -1.7368783  0.9638139
#> anansi_ID_sample_2_1   0.0404813 -1.6560372
#> anansi_ID_sample_3_1  -0.1241795  1.0708611
#> anansi_ID_sample_4_1  -0.6127534 -0.1090264
#> anansi_ID_sample_5_1   0.1606625  1.8991864
#> anansi_ID_sample_6_1  -0.6623596 -1.1370307
#> anansi_ID_sample_7_1  -0.3348517 -0.2797198
#> anansi_ID_sample_8_1   0.6230116 -0.8941291
#> anansi_ID_sample_9_1   1.0281943  0.1367018
#> anansi_ID_sample_10_1 -1.1345782 -0.7491654
#> 
#> [[24]]
#>                               y_6        x_4
#> anansi_ID_sample_1_1   0.39680092  0.9638139
#> anansi_ID_sample_2_1   0.32642296 -1.6560372
#> anansi_ID_sample_3_1  -1.45934331  1.0708611
#> anansi_ID_sample_4_1  -0.79929553 -0.1090264
#> anansi_ID_sample_5_1   0.47001280  1.8991864
#> anansi_ID_sample_6_1  -1.15018052 -1.1370307
#> anansi_ID_sample_7_1  -0.27159672 -0.2797198
#> anansi_ID_sample_8_1   0.45742422 -0.8941291
#> anansi_ID_sample_9_1  -0.01695794  0.1367018
#> anansi_ID_sample_10_1 -0.54159436 -0.7491654
#> 
#> [[25]]
#>                               y_9        x_4
#> anansi_ID_sample_1_1  -1.69018029  0.9638139
#> anansi_ID_sample_2_1  -0.94549766 -1.6560372
#> anansi_ID_sample_3_1  -0.34301196  1.0708611
#> anansi_ID_sample_4_1   0.35846241 -0.1090264
#> anansi_ID_sample_5_1   0.04818886  1.8991864
#> anansi_ID_sample_6_1   1.13017228 -1.1370307
#> anansi_ID_sample_7_1  -0.52763231 -0.2797198
#> anansi_ID_sample_8_1   0.47957533 -0.8941291
#> anansi_ID_sample_9_1   1.41752930  0.1367018
#> anansi_ID_sample_10_1 -0.08009077 -0.7491654
#> 
#> [[26]]
#>                              y_12        x_4
#> anansi_ID_sample_1_1  -1.16985274  0.9638139
#> anansi_ID_sample_2_1   1.32208831 -1.6560372
#> anansi_ID_sample_3_1   1.36596464  1.0708611
#> anansi_ID_sample_4_1  -0.24972371 -0.1090264
#> anansi_ID_sample_5_1  -0.42472135  1.8991864
#> anansi_ID_sample_6_1   1.09843934 -1.1370307
#> anansi_ID_sample_7_1   0.67092304 -0.2797198
#> anansi_ID_sample_8_1   0.03095466 -0.8941291
#> anansi_ID_sample_9_1   2.36851196  0.1367018
#> anansi_ID_sample_10_1  2.58578811 -0.7491654
#> 
#> [[27]]
#>                               y_6         x_5
#> anansi_ID_sample_1_1   0.39680092  0.51819908
#> anansi_ID_sample_2_1   0.32642296 -0.19233721
#> anansi_ID_sample_3_1  -1.45934331  0.02880981
#> anansi_ID_sample_4_1  -0.79929553  0.35859089
#> anansi_ID_sample_5_1   0.47001280 -0.02899503
#> anansi_ID_sample_6_1  -1.15018052  1.14704207
#> anansi_ID_sample_7_1  -0.27159672  0.37358894
#> anansi_ID_sample_8_1   0.45742422  0.32333921
#> anansi_ID_sample_9_1  -0.01695794 -0.82981932
#> anansi_ID_sample_10_1 -0.54159436  1.39446258
#> 
#> [[28]]
#>                               y_7         x_5
#> anansi_ID_sample_1_1   0.87581342  0.51819908
#> anansi_ID_sample_2_1   0.74083815 -0.19233721
#> anansi_ID_sample_3_1  -0.16481239  0.02880981
#> anansi_ID_sample_4_1  -0.75096428  0.35859089
#> anansi_ID_sample_5_1  -1.25336629 -0.02899503
#> anansi_ID_sample_6_1  -1.03623626  1.14704207
#> anansi_ID_sample_7_1  -0.02943857  0.37358894
#> anansi_ID_sample_8_1  -0.27494378  0.32333921
#> anansi_ID_sample_9_1   0.48545147 -0.82981932
#> anansi_ID_sample_10_1  0.91698203  1.39446258
#> 
#> [[29]]
#>                               y_9         x_5
#> anansi_ID_sample_1_1  -1.69018029  0.51819908
#> anansi_ID_sample_2_1  -0.94549766 -0.19233721
#> anansi_ID_sample_3_1  -0.34301196  0.02880981
#> anansi_ID_sample_4_1   0.35846241  0.35859089
#> anansi_ID_sample_5_1   0.04818886 -0.02899503
#> anansi_ID_sample_6_1   1.13017228  1.14704207
#> anansi_ID_sample_7_1  -0.52763231  0.37358894
#> anansi_ID_sample_8_1   0.47957533  0.32333921
#> anansi_ID_sample_9_1   1.41752930 -0.82981932
#> anansi_ID_sample_10_1 -0.08009077  1.39446258
#> 
#> [[30]]
#>                             y_10         x_5
#> anansi_ID_sample_1_1  -0.3941870  0.51819908
#> anansi_ID_sample_2_1   0.5252576 -0.19233721
#> anansi_ID_sample_3_1  -0.7959351  0.02880981
#> anansi_ID_sample_4_1   0.2450964  0.35859089
#> anansi_ID_sample_5_1  -0.7793198 -0.02899503
#> anansi_ID_sample_6_1   0.5949239  1.14704207
#> anansi_ID_sample_7_1   1.1088909  0.37358894
#> anansi_ID_sample_8_1  -0.9429690  0.32333921
#> anansi_ID_sample_9_1   0.7003440 -0.82981932
#> anansi_ID_sample_10_1 -0.4246748  1.39446258
#> 
#> [[31]]
#>                              y_12         x_5
#> anansi_ID_sample_1_1  -1.16985274  0.51819908
#> anansi_ID_sample_2_1   1.32208831 -0.19233721
#> anansi_ID_sample_3_1   1.36596464  0.02880981
#> anansi_ID_sample_4_1  -0.24972371  0.35859089
#> anansi_ID_sample_5_1  -0.42472135 -0.02899503
#> anansi_ID_sample_6_1   1.09843934  1.14704207
#> anansi_ID_sample_7_1   0.67092304  0.37358894
#> anansi_ID_sample_8_1   0.03095466  0.32333921
#> anansi_ID_sample_9_1   2.36851196 -0.82981932
#> anansi_ID_sample_10_1  2.58578811  1.39446258
#> 
#> [[32]]
#>                              y_1         x_6
#> anansi_ID_sample_1_1   0.8250203 -0.19154358
#> anansi_ID_sample_2_1  -0.7501828  0.27227702
#> anansi_ID_sample_3_1   0.5362105 -1.08165901
#> anansi_ID_sample_4_1  -1.5727175 -2.32939936
#> anansi_ID_sample_5_1  -0.9864274 -0.54962596
#> anansi_ID_sample_6_1   1.9835123 -0.07257999
#> anansi_ID_sample_7_1  -1.8513784  1.03228840
#> anansi_ID_sample_8_1  -0.9097882  0.21513853
#> anansi_ID_sample_9_1  -1.9514400 -0.49434012
#> anansi_ID_sample_10_1 -0.8002765  1.51213828
#> 
#> [[33]]
#>                               y_2         x_6
#> anansi_ID_sample_1_1  -1.86936160 -0.19154358
#> anansi_ID_sample_2_1  -0.75057148  0.27227702
#> anansi_ID_sample_3_1  -0.59093839 -1.08165901
#> anansi_ID_sample_4_1  -0.74209927 -2.32939936
#> anansi_ID_sample_5_1   0.69351208 -0.54962596
#> anansi_ID_sample_6_1  -0.05946397 -0.07257999
#> anansi_ID_sample_7_1  -1.86389510  1.03228840
#> anansi_ID_sample_8_1  -1.27450892  0.21513853
#> anansi_ID_sample_9_1  -1.78177081 -0.49434012
#> anansi_ID_sample_10_1 -0.50859339  1.51213828
#> 
#> [[34]]
#>                              y_4         x_6
#> anansi_ID_sample_1_1   0.9169111 -0.19154358
#> anansi_ID_sample_2_1   1.2138748  0.27227702
#> anansi_ID_sample_3_1  -0.6859325 -1.08165901
#> anansi_ID_sample_4_1  -1.0945682 -2.32939936
#> anansi_ID_sample_5_1   0.3648741 -0.54962596
#> anansi_ID_sample_6_1  -1.0080171 -0.07257999
#> anansi_ID_sample_7_1   0.5549964  1.03228840
#> anansi_ID_sample_8_1   0.4360659  0.21513853
#> anansi_ID_sample_9_1   0.1053733 -0.49434012
#> anansi_ID_sample_10_1 -0.2537379  1.51213828
#> 
#> [[35]]
#>                              y_5         x_6
#> anansi_ID_sample_1_1  -0.9764912 -0.19154358
#> anansi_ID_sample_2_1   0.3942178  0.27227702
#> anansi_ID_sample_3_1   1.0149030 -1.08165901
#> anansi_ID_sample_4_1   0.6314438 -2.32939936
#> anansi_ID_sample_5_1  -0.4833983 -0.54962596
#> anansi_ID_sample_6_1  -1.5263339 -0.07257999
#> anansi_ID_sample_7_1   0.6168423  1.03228840
#> anansi_ID_sample_8_1   1.2292868  0.21513853
#> anansi_ID_sample_9_1   0.1574530 -0.49434012
#> anansi_ID_sample_10_1  1.4063421  1.51213828
#> 
#> [[36]]
#>                               y_9         x_6
#> anansi_ID_sample_1_1  -1.69018029 -0.19154358
#> anansi_ID_sample_2_1  -0.94549766  0.27227702
#> anansi_ID_sample_3_1  -0.34301196 -1.08165901
#> anansi_ID_sample_4_1   0.35846241 -2.32939936
#> anansi_ID_sample_5_1   0.04818886 -0.54962596
#> anansi_ID_sample_6_1   1.13017228 -0.07257999
#> anansi_ID_sample_7_1  -0.52763231  1.03228840
#> anansi_ID_sample_8_1   0.47957533  0.21513853
#> anansi_ID_sample_9_1   1.41752930 -0.49434012
#> anansi_ID_sample_10_1 -0.08009077  1.51213828
#> 
#> [[37]]
#>                             y_11         x_6
#> anansi_ID_sample_1_1  -1.1431376 -0.19154358
#> anansi_ID_sample_2_1   0.2334402  0.27227702
#> anansi_ID_sample_3_1  -0.1892366 -1.08165901
#> anansi_ID_sample_4_1  -1.6634113 -2.32939936
#> anansi_ID_sample_5_1   1.9156908 -0.54962596
#> anansi_ID_sample_6_1  -0.8168944 -0.07257999
#> anansi_ID_sample_7_1   0.3836505  1.03228840
#> anansi_ID_sample_8_1  -0.4586323  0.21513853
#> anansi_ID_sample_9_1  -0.7157040 -0.49434012
#> anansi_ID_sample_10_1  0.4295485  1.51213828
#> 
#> [[38]]
#>                              y_12         x_6
#> anansi_ID_sample_1_1  -1.16985274 -0.19154358
#> anansi_ID_sample_2_1   1.32208831  0.27227702
#> anansi_ID_sample_3_1   1.36596464 -1.08165901
#> anansi_ID_sample_4_1  -0.24972371 -2.32939936
#> anansi_ID_sample_5_1  -0.42472135 -0.54962596
#> anansi_ID_sample_6_1   1.09843934 -0.07257999
#> anansi_ID_sample_7_1   0.67092304  1.03228840
#> anansi_ID_sample_8_1   0.03095466  0.21513853
#> anansi_ID_sample_9_1   2.36851196 -0.49434012
#> anansi_ID_sample_10_1  2.58578811  1.51213828
#> 
#> [[39]]
#>                              y_1        x_7
#> anansi_ID_sample_1_1   0.8250203 -0.5794633
#> anansi_ID_sample_2_1  -0.7501828  1.6747896
#> anansi_ID_sample_3_1   0.5362105 -1.0009803
#> anansi_ID_sample_4_1  -1.5727175  1.2227028
#> anansi_ID_sample_5_1  -0.9864274  1.0771882
#> anansi_ID_sample_6_1   1.9835123 -0.6119455
#> anansi_ID_sample_7_1  -1.8513784  0.5068670
#> anansi_ID_sample_8_1  -0.9097882  0.4600684
#> anansi_ID_sample_9_1  -1.9514400  1.4843919
#> anansi_ID_sample_10_1 -0.8002765  0.8819632
#> 
#> [[40]]
#>                               y_2        x_7
#> anansi_ID_sample_1_1  -1.86936160 -0.5794633
#> anansi_ID_sample_2_1  -0.75057148  1.6747896
#> anansi_ID_sample_3_1  -0.59093839 -1.0009803
#> anansi_ID_sample_4_1  -0.74209927  1.2227028
#> anansi_ID_sample_5_1   0.69351208  1.0771882
#> anansi_ID_sample_6_1  -0.05946397 -0.6119455
#> anansi_ID_sample_7_1  -1.86389510  0.5068670
#> anansi_ID_sample_8_1  -1.27450892  0.4600684
#> anansi_ID_sample_9_1  -1.78177081  1.4843919
#> anansi_ID_sample_10_1 -0.50859339  0.8819632
#> 
#> [[41]]
#>                               y_6        x_7
#> anansi_ID_sample_1_1   0.39680092 -0.5794633
#> anansi_ID_sample_2_1   0.32642296  1.6747896
#> anansi_ID_sample_3_1  -1.45934331 -1.0009803
#> anansi_ID_sample_4_1  -0.79929553  1.2227028
#> anansi_ID_sample_5_1   0.47001280  1.0771882
#> anansi_ID_sample_6_1  -1.15018052 -0.6119455
#> anansi_ID_sample_7_1  -0.27159672  0.5068670
#> anansi_ID_sample_8_1   0.45742422  0.4600684
#> anansi_ID_sample_9_1  -0.01695794  1.4843919
#> anansi_ID_sample_10_1 -0.54159436  0.8819632
#> 
#> [[42]]
#>                               y_7        x_7
#> anansi_ID_sample_1_1   0.87581342 -0.5794633
#> anansi_ID_sample_2_1   0.74083815  1.6747896
#> anansi_ID_sample_3_1  -0.16481239 -1.0009803
#> anansi_ID_sample_4_1  -0.75096428  1.2227028
#> anansi_ID_sample_5_1  -1.25336629  1.0771882
#> anansi_ID_sample_6_1  -1.03623626 -0.6119455
#> anansi_ID_sample_7_1  -0.02943857  0.5068670
#> anansi_ID_sample_8_1  -0.27494378  0.4600684
#> anansi_ID_sample_9_1   0.48545147  1.4843919
#> anansi_ID_sample_10_1  0.91698203  0.8819632
#> 
#> [[43]]
#>                             y_11        x_7
#> anansi_ID_sample_1_1  -1.1431376 -0.5794633
#> anansi_ID_sample_2_1   0.2334402  1.6747896
#> anansi_ID_sample_3_1  -0.1892366 -1.0009803
#> anansi_ID_sample_4_1  -1.6634113  1.2227028
#> anansi_ID_sample_5_1   1.9156908  1.0771882
#> anansi_ID_sample_6_1  -0.8168944 -0.6119455
#> anansi_ID_sample_7_1   0.3836505  0.5068670
#> anansi_ID_sample_8_1  -0.4586323  0.4600684
#> anansi_ID_sample_9_1  -0.7157040  1.4843919
#> anansi_ID_sample_10_1  0.4295485  0.8819632
#> 
#> [[44]]
#>                              y_12        x_7
#> anansi_ID_sample_1_1  -1.16985274 -0.5794633
#> anansi_ID_sample_2_1   1.32208831  1.6747896
#> anansi_ID_sample_3_1   1.36596464 -1.0009803
#> anansi_ID_sample_4_1  -0.24972371  1.2227028
#> anansi_ID_sample_5_1  -0.42472135  1.0771882
#> anansi_ID_sample_6_1   1.09843934 -0.6119455
#> anansi_ID_sample_7_1   0.67092304  0.5068670
#> anansi_ID_sample_8_1   0.03095466  0.4600684
#> anansi_ID_sample_9_1   2.36851196  1.4843919
#> anansi_ID_sample_10_1  2.58578811  0.8819632
#> 
#> [[45]]
#>                              y_1        x_8
#> anansi_ID_sample_1_1   0.8250203 -0.5367022
#> anansi_ID_sample_2_1  -0.7501828  1.2855537
#> anansi_ID_sample_3_1   0.5362105  0.5878485
#> anansi_ID_sample_4_1  -1.5727175 -1.3084820
#> anansi_ID_sample_5_1  -0.9864274  0.3167263
#> anansi_ID_sample_6_1   1.9835123  1.1941583
#> anansi_ID_sample_7_1  -1.8513784  0.9130571
#> anansi_ID_sample_8_1  -0.9097882 -0.7867530
#> anansi_ID_sample_9_1  -1.9514400 -0.4109297
#> anansi_ID_sample_10_1 -0.8002765  0.4763737
#> 
#> [[46]]
#>                              y_3        x_8
#> anansi_ID_sample_1_1  -1.7368783 -0.5367022
#> anansi_ID_sample_2_1   0.0404813  1.2855537
#> anansi_ID_sample_3_1  -0.1241795  0.5878485
#> anansi_ID_sample_4_1  -0.6127534 -1.3084820
#> anansi_ID_sample_5_1   0.1606625  0.3167263
#> anansi_ID_sample_6_1  -0.6623596  1.1941583
#> anansi_ID_sample_7_1  -0.3348517  0.9130571
#> anansi_ID_sample_8_1   0.6230116 -0.7867530
#> anansi_ID_sample_9_1   1.0281943 -0.4109297
#> anansi_ID_sample_10_1 -1.1345782  0.4763737
#> 
#> [[47]]
#>                              y_4        x_8
#> anansi_ID_sample_1_1   0.9169111 -0.5367022
#> anansi_ID_sample_2_1   1.2138748  1.2855537
#> anansi_ID_sample_3_1  -0.6859325  0.5878485
#> anansi_ID_sample_4_1  -1.0945682 -1.3084820
#> anansi_ID_sample_5_1   0.3648741  0.3167263
#> anansi_ID_sample_6_1  -1.0080171  1.1941583
#> anansi_ID_sample_7_1   0.5549964  0.9130571
#> anansi_ID_sample_8_1   0.4360659 -0.7867530
#> anansi_ID_sample_9_1   0.1053733 -0.4109297
#> anansi_ID_sample_10_1 -0.2537379  0.4763737
#> 
#> [[48]]
#>                               y_9        x_8
#> anansi_ID_sample_1_1  -1.69018029 -0.5367022
#> anansi_ID_sample_2_1  -0.94549766  1.2855537
#> anansi_ID_sample_3_1  -0.34301196  0.5878485
#> anansi_ID_sample_4_1   0.35846241 -1.3084820
#> anansi_ID_sample_5_1   0.04818886  0.3167263
#> anansi_ID_sample_6_1   1.13017228  1.1941583
#> anansi_ID_sample_7_1  -0.52763231  0.9130571
#> anansi_ID_sample_8_1   0.47957533 -0.7867530
#> anansi_ID_sample_9_1   1.41752930 -0.4109297
#> anansi_ID_sample_10_1 -0.08009077  0.4763737
#> 

pairwiseApply(
    X = web,
    FUN = function(x, y) cor(x, y),
    MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE
)
#>       x_1y_1       x_1y_3       x_1y_4       x_1y_9      x_1y_11      x_1y_12 
#>  0.522113131 -0.696267873 -0.531100148 -0.009102863 -0.459253107 -0.248220130 
#>       x_2y_1       x_2y_2       x_2y_4       x_2y_5       x_2y_8       x_2y_9 
#>  0.352041701  0.486968863 -0.499887439 -0.002513495  0.428928257  0.192281113 
#>      x_2y_11       x_3y_2       x_3y_4       x_3y_6       x_3y_7       x_3y_8 
#> -0.096951314 -0.017878495  0.161562745 -0.064485375  0.026846110 -0.063464116 
#>      x_3y_10      x_3y_12       x_4y_1       x_4y_2       x_4y_3       x_4y_6 
#>  0.573111776  0.296869501 -0.007428255  0.168376016 -0.052740480  0.062542941 
#>       x_4y_9      x_4y_12       x_5y_6       x_5y_7       x_5y_9      x_5y_10 
#> -0.182695449 -0.420624713 -0.348240559 -0.040992676 -0.105253848 -0.129788469 
#>      x_5y_12       x_6y_1       x_6y_2       x_6y_4       x_6y_5       x_6y_9 
#> -0.012472361 -0.012763393 -0.175092409  0.474790609  0.157377526 -0.199612411 
#>      x_6y_11      x_6y_12       x_7y_1       x_7y_2       x_7y_6       x_7y_7 
#>  0.444308301  0.386922730 -0.797650050  0.020269444  0.484813682  0.133899906 
#>      x_7y_11      x_7y_12       x_8y_1       x_8y_3       x_8y_4       x_8y_9 
#>  0.246368329  0.227181472  0.335983986 -0.047652392  0.144181990 -0.143731386 
```
