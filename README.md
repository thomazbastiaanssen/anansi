<!-- README.md is generated from README.Rmd. Please edit that file -->

The `anansi` package computes and compares the association between the
features of two ’omics datasets that are known to interact based on a
database such as KEGG.

If you use this software, please cite our work.

``` r
citation("anansi")
```

    ## Warning in citation("anansi"): no date field in DESCRIPTION file of package
    ## 'anansi'

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

OK, now let’s get started.

``` r
#install and load anansi
#devtools::install_github("thomazbastiaanssen/anansi")
library(anansi)

#load ggplot2 and ggforce to plot results
library(ggplot)
library(ggforce)

#load anansi dictionary and human-readable names for KEGG compounds and orthologues
data(dictionary)
data(cpd_translation)
data(KO_translation)

#load example data from FMT Aging study
data(FMT_metadata)
data(FMT_metabs)
data(FMT_KOs)
```

## Data preparation

The main anansi function expects data in the `anansiWeb` format;
Basically a list with exactly three tables: The first table, `tableY`,
should be a count table of metabolites. The second table, `tableX`,
should be a count table of functions. Both tables should have rows as
features and columns as samples.

the third table should be a binary adjacency matrix with the column
names of table Y as rows and the colnames of table X as columns. Such a
table is provided in the anansi library and is referred to as a
dictionary (because you use it to look up which metabolites interact
with which functions).

the `weaveWebFromTables()` function can be used to parse these tables
into an `anansiWeb` object.

``` r
#Clean and CLR-transform the KEGG orthologue table
KOs   <- floor(KOs)
KOs   <- apply(KOs,c(1,2),function(x) as.numeric(as.character(x)))
KOs   <- KOs[apply(KOs == 0, 1, sum) <= (ncol(KOs) * 0.90), ] 

#only keep functions that are represented in the dictionary
KOs   <- KOs[row.names(KOs) %in% sort(unique(unlist(anansi_dic))),]

KOs.exp = clr_lite(KOs)

#anansi expects samples to be rows and features to be columns. 
t1 = t(metab)
t2 = t(KOs.exp)
```

## weave a web

``` r
web = weaveWebFromTables(tableY = t1, tableX = t2, dictionary = anansi_dic)
```

## Run anansi

``` r
anansi_out = anansi(web    = web, #generated above
                    method = "pearson", #define the type of correlation used
                    groups = metadata$Legend, #optional, to compare associations between groups
                    adjust.method = "BH", #apply the Benjamini-Hochberg procedure for FDR
                    verbose = T #To let you know what's happening
                    )
```

Anansi gives a complex nested list of lists as an output. Two functions
exist that will wrangle your data to more friendly formats for you. You
can either use spinToLong() or spinToWide(). They will give you long or
wide format data.frames, respectively. For general reporting, we
recommend sticking to wide format. Long format can be helpful to plug
the data into ggplot2.

``` r
anansiLong = spinToLong(anansi_output = anansi_out)  

#Now it's ready to be plugged into ggplot2, though let's clean up a bit more. 

#Only consider interactions where the entire model fits well enough
anansiLong = anansiLong[anansiLong$model_full_q.values < 0.1,]

#Plug into ggplot2

ggplot(data = anansiLong, 
       aes(x      = estimates, 
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
  scale_fill_manual(values = c("correlation_Young yFMT" = "#2166ac", 
                               "correlation_Aged oFMT"  = "#b2182b", 
                               "correlation_Aged yFMT"  = "#ef8a62", 
                               "correlation_All"        = "gray"))+
  theme_bw() + 
  ylab("") + 
  xlab("Pearson's rho")
```
