
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Knowledge-based multi-modal integration using anansi <img src="man/figures/anansi_hex.png" align="right" width="120" alt="The anansi hex sticker" />

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/thomazbastiaanssen/anansi)](https://github.com/thomazbastiaanssen/anansi/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/thomazbastiaanssen/anansi)](https://github.com/thomazbastiaanssen/anansi/pulls)
<!-- badges: end -->

## Introduction

The `anansi` package computes and compares the association between the
features of two ’omics data sets that are known to interact based on a
database such as KEGG. Studies including both microbiome and
metabolomics data are becoming more common. Often, it would be helpful
to integrate both data sets in order to see if they corroborate each
others patterns. All vs all association is imprecise and likely to yield
spurious associations. This package takes a knowledge-based approach to
constrain association search space, only considering metabolite-function
interactions that have been recorded in a pathway database. This package
also provides a framework to assess differential association.

While `anansi` is geared towards metabolite-function interactions in the
context of host-microbe interactions, it is perfectly capable of
handling any other pair of data sets where some features interact
canonically.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `anansi` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("anansi")
```

And the development version from
[GitHub](https://github.com/thomazbastiaanssen/anansi) with `remotes`:

``` r
install.packages("remotes")
remotes::install_github("thomazbastiaanssen/anansi")
```

## Getting started using anansi

[See the vignettes on the package
site.](https://thomazbastiaanssen.github.io/anansi/articles/anansi.html)
Additionally, a [preprint](https://arxiv.org/abs/2305.10832) is
available.

## Citation

Below is the citation output from using `citation('anansi')` in R.

``` r
print(citation('anansi'), bibtex = TRUE)
#> To cite anansi in publications use:
#> 
#>   Bastiaanssen TFS, Quinn TP, Cryan JF (2023) Knowledge-based
#>   Integration of Multi-Omic Datasets with Anansi: Annotation-based
#>   Analysis of Specific Interactions arXiv. doi:
#>   10.48550/arXiv.2305.10832
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Knowledge-based Integration of Multi-Omic Datasets with Anansi: Annotation-based Analysis of Specific Interactions},
#>     author = {Thomaz F S Bastiaanssen and Thomas P Quinn and John F Cryan},
#>     journal = {arXiv},
#>     year = {2023},
#>     doi = {10.48550/arXiv.2305.10832},
#>   }
```

## Code of Conduct

Please note that the `anansi` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.
