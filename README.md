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
    ##   person and comment = c) (2022). anansi: ANnotation based ANalysis of
    ##   Specific Interactions. R package version 0.5.0.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {anansi: ANnotation based ANalysis of Specific Interactions},
    ##     author = {{person} and comment = c)},
    ##     year = {2022},
    ##     note = {R package version 0.5.0},
    ##   }
    ## 
    ## ATTENTION: This citation information has been auto-generated from the
    ## package DESCRIPTION file and may need manual editing, see
    ## 'help("citation")'.

OK, now let’s get started.

``` r
#install and load anansi
#devtools::install_github("thomazbastiaanssen/anansi")
library(anansi)

#load tidyverse, whighly recommended for data wrangling 
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.6     ✓ dplyr   1.0.7
    ## ✓ tidyr   1.1.4     ✓ stringr 1.4.0
    ## ✓ readr   2.0.0     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
#load ggforce for plotting results
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

The main anansi function expects data in the ‘web’ format; a list with
exactly three elements: The first table should be a count table of
metabolites. The second table should be a count table of functions. Both
tables should have rows as features and columns as samples.

the third table should be a binary adjacency matrix with the colnames of
table 1 as rows and the colnames of table 2 as columns. Such a table is
provided in the anansi library and is referred to as a dictionary
(becasue you use it to look up which metabolites interact with which
functions).

the weaveWebFromTables() function can be used to parse these tables into

``` r
#First rid of the "cpd:" part of the names of compounds to match the metabolite table
names(anansi_dic) = gsub(pattern = "cpd:", replacement = "", x = names(anansi_dic))

#Then prune the dictionary list to reduce memory usage
anansi_dic = anansi_dic[names(anansi_dic) %in% rownames(metab)]

#then wrangle the adjacency list to a binary matrix
anansi_dic <- anansi_dic %>% 
  enframe() %>% 
  unnest(cols = c(value)) %>% add_column(present = T) %>%
  pivot_wider(values_from = present, values_fill = F) %>%
  filter(!is.na(value)) %>% 
  column_to_rownames("value") %>% 
  t()

#Clean and CLR-transform the KEGG orthologue table
KOs   <- floor(KOs)
KOs   <- apply(KOs,c(1,2),function(x) as.numeric(as.character(x)))
KOs   <- KOs[apply(KOs == 0, 1, sum) <= (ncol(KOs) * 0.90), ] 

#only keep functions that are represented in the dictionary
KOs   <- KOs[row.names(KOs) %in% colnames(anansi_dic),]

KOs.exp = clr_lite(KOs)

#anansi expects samples to be rows and features to be columns. 
t1 = t(metab)
t2 = t(KOs.exp)
```

## weave a web

``` r
web = weaveWebFromTables(table1 = t1, table2 = t2, dictionary = anansi_dic)
```

    ## [1] "3 were matched between table 1 and the columns of the adjacency matrix"
    ## [1] "50 were matched between table 2 and the rows of the adjacency matrix"

## Run anansi

``` r
anansi_out = anansi(web    = web, #generated above
                    method = "pearson", #define the type of correlation used
                    groups = metadata$Legend, #optional, to compare associations between groups
                    adjust.method = "BH", #apply the Benjamini-Hochberg procedure for FDR
                    verbose = T #To let you know what's happening
                    )
```

    ## [1] "Running annotation-based correlations"
    ## [1] "Running correlations for the following groups: Aged yFMT, Aged oFMT, Young yFMT and all together"
    ## [1] "Fitting models for differential correlation testing"

Anansi gives a nested list of lists as an output. I highly recommend
using purrr from the tidyverse to parse it

``` r
anansi_out %>%
  
  #Induvidually wrangle all results to long format
  map_depth(3, ~data.frame(.x)) %>%
  map_depth(3, ~rownames_to_column(.x, "Compound")) %>%
  map_depth(3, ~pivot_longer(.x, !Compound, names_to = "Function")) %>%
  flatten()  %>%
  
  #Merge estimates, p and q values per model
  map(~reduce(.x, full_join, by = c("Compound", "Function"))) %>%
  map(~rename(.x, c("Estimate" = "value.x", "p.value" = "value.y", "q.value" = "value"))) %>%
  
  #Merge all data to one data.frame
  bind_rows(.id = "Type") %>%
  
  #Organise data by model and parameter type
  pivot_wider(names_from = Type, values_from = c(Estimate, p.value, q.value)) %>%
  
  #Handle correlations separately from total model parameters
  pivot_longer(cols = -c(Compound, Function, 
                         Estimate_modelfit,   p.value_modelfit,   q.value_modelfit, 
                         Estimate_emergent,   p.value_emergent,   q.value_emergent, 
                         Estimate_disjointed, p.value_disjointed, q.value_disjointed)) %>% 
  separate(col = name, into = c("Parameter", "Group"), sep = "_") %>% 
  pivot_wider(names_from = Parameter) %>% 
  #Now it's ready to be piped into ggplot2, though let's clean up a bit more. 
  
  #Only consider canonical interactions
  filter(p.value != 1) %>%
  
  #Only consider interactions where the entire model fits well enough
  filter(q.value_modelfit < 0.1) %>% 
  
  group_by(Compound, Function) %>% 
  filter(min(p.value) < 0.05) %>%
  ungroup() %>%
  
  #Make names of KEGG compounds and functions human-readable
  inner_join(KO_translation, by = "Function") %>%
  inner_join(cpd_translation, by = "Compound") %>%
  
  #Make names pretty
  mutate(Function = paste(Function, KO_name, sep = ": "), 
         Compound = paste(Compound, Compound_name, sep = ": "), 
         Compound = str_remove(Compound, pattern = ";.*")) %>% 
  
  
  #Pipe into ggplot2
  ggplot(aes(x      = Estimate, 
             y      = Function, 
             fill   = Group, 
             alpha  = p.value_disjointed < 0.05)) + 
  
  #Make a vertical dashed red line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red")+
  
  #Points show  raw correlation coefficients
  geom_point(shape = 21, size = 3) + 
  
  #facet per compound
  ggforce::facet_col(~Compound, space = "free", scales = "free_y") + 
  
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

![](README_files/figure-markdown_github/unnamed-chunk-6-1.png)
