<!-- README.md is generated from README.Rmd. Please edit that file -->

## Introduction

The `anansi` package computes and compares the association between the
features of two ’omics datasets that are known to interact based on a
database such as KEGG. Studies including both microbiome and
metabolomics data are becoming more common. Often, it would be helpful
to integrate both datasets in order to see if they corroborate each
others patterns. All vs all association is imprecise and likely to yield
spurious associations. This package takes a knowledge-based approach to
constrain association search space, only considering metabolite-function
interactions that have been recorded in a pathway database. This package
also provides a framework to assess differential association.

While `anansi` is geared towards metabolite-function interactions in the
context of host-microbe interactions, it is perfectly capable of
handling any other pair of datasets where some features interact
canonically.

If you use this software, please cite our work.

``` r
citation("anansi")
```

    ## 
    ## To cite package 'anansi' in publications use:
    ## 
    ##   Thomaz Bastiaanssen (2022). anansi: Annotation-based Analysis of
    ##   Specific Interactions. R package version 0.5.0.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {anansi: Annotation-based Analysis of Specific Interactions},
    ##     author = {Thomaz Bastiaanssen},
    ##     year = {2022},
    ##     note = {R package version 0.5.0},
    ##   }

## Setup

OK, now let’s get started. We’ll load a complementary training dataset
using `data(FMT_data)`. This loads a curated snippet from the dataset
described in more detail here:
<https://doi.org/10.1038/s43587-021-00093-9>  
A very early version of `anansi` was used to generate “Extended Data
Fig. 7” in that paper.

``` r
#install and load anansi
#devtools::install_github("thomazbastiaanssen/anansi")
library(anansi)

#load ggplot2 and ggforce to plot results
library(ggplot2)
library(ggforce)

#load anansi dictionary and complementary human-readable names for KEGG compounds and orthologues
data(dictionary)

#load example data + metadata from FMT Aging study
data(FMT_data)
```

## Data preparation

The main `anansi` function expects data in the `anansiWeb` format;
Basically a list with exactly three tables: The first table, `tableY`,
should be a count table of metabolites. The second table, `tableX`,
should be a count table of functions. Both tables should have rows as
features and columns as samples.

the third table should be a binary adjacency matrix with the column
names of `tableY` as rows and the column names of `tableX` as columns.
Such an adjacency matrix is provided in the `anansi` library and is
referred to as a dictionary (because you use it to look up which
metabolites interact with which functions).

### A note on functional microbiome data

Two common questions in the host-microbiome field are “Who’s there?” and
“What are they doing?”. Techniques like 16S sequencing and shotgun
metagenomics sequencing are most commonly used to answer the first
question. The second question can be a bit more tricky - often we’ll
need functional inference software to address them. For 16S sequencing,
algorithms like PICRUSt2 and Piphillin can be used to infer function.
For shotgun metagenomics, HUMANn3 in the bioBakery suite can be used.  
All of these algorithms can produce functional count data in terms of
KEGG Orthologues (KOs). These tables can be directly plugged in to
`anansi`.

``` r
#Clean and CLR-transform the KEGG orthologue table.

#Only keep functions that are represented in the dictionary.
KOs     <- FMT_KOs[row.names(FMT_KOs) %in% sort(unique(unlist(anansi_dic))),]

#Cut the decimal part off.
KOs     <- floor(KOs)

#Ensure all entires are numbers.
KOs     <- apply(KOs,c(1,2),function(x) as.numeric(as.character(x)))

#Remove all features with < 10% prevalence in the dataset.
KOs     <- KOs[apply(KOs == 0, 1, sum) <= (ncol(KOs) * 0.90), ] 

#Perform a centered log-ratio transformation on the functional count table.
KOs.exp <- clr_lite(KOs)

#anansi expects samples to be rows and features to be columns. 
t1      <- t(FMT_metab)
t2      <- t(KOs.exp)
```

## Weave a web🕸️

The `weaveWebFromTables()` function can be used to parse the tables that
we prepared above into an `anansiWeb` object. The `anansiWeb` format is
a necessary input file for the main `anansi` workflow.

``` r
web <- weaveWebFromTables(tableY = t1, tableX = t2, dictionary = anansi_dic)
```

    ## [1] "3 were matched between table 1 and the columns of the adjacency matrix"
    ## [1] "50 were matched between table 2 and the rows of the adjacency matrix"

## Run anansi🕷️

The main workspider in this package is called `anansi`. Generally, you
want to give it two arguments. First, there’s `web`, which is an
`ananisWeb` object, such as the one we generated in the above step.
Second, there’s `groups`, which should be a vector to compare the
associations on. For instance, this may be a vector containing
categories such as your treatment groups, or even a continuous value
like age or .

``` r
anansi_out <- anansi(web    = web, #generated above
                     method = "pearson", #define the type of correlation used
                     groups = FMT_metadata$Legend, #Compare associations between treatments
                     adjust.method = "BH", #apply the Benjamini-Hochberg procedure for FDR
                     verbose = T #To let you know what's happening
                     )
```

    ## [1] "Running annotation-based correlations"
    ## [1] "Running correlations for the following groups: Aged yFMT, Aged oFMT, Young yFMT and all together"
    ## [1] "Fitting models for differential correlation testing"

## Spin to a table📝

`anansi` gives a complex nested list of lists as an output. Two
functions exist that will wrangle your data to more friendly formats for
you. You can either use `spinToLong()` or `spinToWide()`. They will give
you long or wide format data.frames, respectively. For general
reporting, we recommend sticking to the wide format as it’s the most
legible.

``` r
anansiLong <- spinToLong(anansi_output = anansi_out)  
#Now it's ready to be plugged into ggplot2, though let's clean up a bit more. 

#Only consider interactions where the entire model fits well enough. 
anansiLong <- anansiLong[anansiLong$model_full_q.values < 0.1,]
```

## Plot the results

The long format can be helpful to plug the data into `ggplot2`. Here, we
recreate part of the results from the FMT Aging study.

``` r
ggplot(data = anansiLong, 
       aes(x      = r.values, 
           y      = Functions, 
           fill   = type, 
           alpha  = model_disjointed_p.values < 0.05)) + 
  
  #Make a vertical dashed red line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red")+
  
  #Points show  raw correlation coefficients
  geom_point(shape = 21, size = 3) + 
  
  #facet per compound
  ggforce::facet_col(~Compounds, space = "free", scales = "free_y") + 
  
  #fix the scales, labels, theme and other layout
  scale_y_discrete(limits = rev, position = "right") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1/3)) +
  scale_fill_manual(values = c("Young yFMT" = "#2166ac", 
                               "Aged oFMT"  = "#b2182b", 
                               "Aged yFMT"  = "#ef8a62", 
                               "All"        = "gray"))+
  theme_bw() + 
  ylab("") + 
  xlab("Pearson's rho")
```

![](README_files/figure-markdown_github/plot_FMT-1.png)
