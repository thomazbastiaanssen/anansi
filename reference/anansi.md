# anansi: Annotation-based Analysis of Specific Interactions.

The anansi package package provides tools to prepare and facilitate
integrative association analysis between the features of two data sets
that are known to interact.

### 1. Input for `anansi()` with [`AnansiWeb()`](https://thomazbastiaanssen.github.io/anansi/reference/AnansiWeb-class.md) and [`MultiFactor()`](https://thomazbastiaanssen.github.io/anansi/reference/MultiFactor-class.md)

- [randomAnansi](https://thomazbastiaanssen.github.io/anansi/reference/randomAnansi.md),
  [`kegg_link()`](https://thomazbastiaanssen.github.io/anansi/reference/kegg_link.md):
  Generate example input

- [AnansiWeb](https://thomazbastiaanssen.github.io/anansi/reference/AnansiWeb-class.md),
  [MultiFactor](https://thomazbastiaanssen.github.io/anansi/reference/MultiFactor-class.md):
  Handle and manipulate input

### 2. Output and cross-compatibility

- [`getGraph()`](https://thomazbastiaanssen.github.io/anansi/reference/getGraph.md):
  Compatibility with igraph::igraph

- [`plotAnansi()`](https://thomazbastiaanssen.github.io/anansi/reference/plotAnansi.md):
  Plot output in the style of miaViz::miaViz

### 3. Vignettes

- [1. Getting started with
  anansi](https://thomazbastiaanssen.github.io/anansi/articles/anansi.html)

- [2. Adjacency
  matrices](https://thomazbastiaanssen.github.io/anansi/articles/adjacency_matrices.html)

- [3. Association
  testing](https://thomazbastiaanssen.github.io/anansi/articles/association_testing.html)

### Description: The `anansi()` function

The main workspider function in the anansi package is called `anansi`.
It manages the individual functionalities of anansi, including
correlation analysis, correlation by group and differential
associations.

## Usage

``` r
anansi(
  web,
  formula,
  groups = NULL,
  metadata = NULL,
  adjust.method = "BH",
  verbose = TRUE,
  return.format = "table",
  ...
)
```

## Arguments

- web:

  An `AnansiWeb` object, containing two tables with 'omics data and a
  dictionary that links them. See `weaveWebFromTables()` for how to
  weave a web.

- formula:

  A formula object. Used to assess differential associations.

- groups:

  A vector of the column names of categorical values in the metadata
  object to control which groups should be assessed for simple
  correlations. If no argument provided, anansi will let you know and
  still to regular correlations according to your dictionary.

- metadata:

  A vector or data.frame of categorical or continuous value necessary
  for differential correlations. Typically a state or treatment. If no
  argument provided, anansi will let you know and still to regular
  correlations according to your dictionary.

- adjust.method:

  Method to adjust p-values for multiple comparisons.
  `adjust.method = "BH"` is the default value. See
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) in the base R
  `stats` package.

- verbose:

  `Logical scalar`. Whether to print diagnostic information (Default:
  `TRUE`). while running. Useful for debugging errors on large datasets.

- return.format:

  `Character scalar`. Should be one of `"table"`, `"list"`, or `"raw"`.
  Should the output of `anansi()` respectively be a wide `data.frame` of
  results, a list containing the results and input, or a list of raw
  output (used for testing purposes). convenient use. (Default:
  `"table"`)

- ...:

  additional arguments (currently not used).

## Value

A list of lists containing correlation coefficients, p-values and
q-values for all operations.

## See also

Useful links:

- <https://github.com/thomazbastiaanssen/anansi>

- <https://thomazbastiaanssen.github.io/anansi>

- Report bugs at <https://github.com/thomazbastiaanssen/anansi/issues>

## Author

**Maintainer**: Thomaz Bastiaanssen <thomazbastiaanssen@gmail.com>
([ORCID](https://orcid.org/0000-0001-6891-734X))

Authors:

- Thomas Quinn <contacttomquinn@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-0286-6329))

- Giulio Benedetti <giulio.benedetti@utu.fi>
  ([ORCID](https://orcid.org/0000-0002-8732-7692))

- Tuomas Borman <tuomas.v.borman@utu.fi>
  ([ORCID](https://orcid.org/0000-0002-8563-8884))

- Leo Lahti <leo.lahti@utu.fi>
  ([ORCID](https://orcid.org/0000-0001-5537-637X))

## Examples

``` r
# Load example data:

data(FMT_data)

# Make sure that columns are features and rows are samples.

t1 <- t(FMT_metab)
t2 <- t(FMT_KOs)

# Run anansi pipeline.

web <- weaveWeb(
    cpd ~ ko,
    tableY = t1,
    tableX = t2,
    link = kegg_link()
)
#> Dropped features in tableX: 2476 remain. 
#> Warning: Argument `metadata` not provided; Please validate sample ID order.

anansi_out <- anansi(
    web = web,
    formula = ~Legend,
    groups = "Legend",
    metadata = FMT_metadata,
    adjust.method = "BH",
    verbose = TRUE
)
#> Fitting least-squares for following model:
#> ~ x + Legend + x:Legend 
#> Running correlations for the following groups:
#>  Aged yFMT, Aged oFMT, Young yFMT


library(tidyr)

# Use tidyr to wrangle the correlation r-values to a single column
anansiLong <- anansi_out |>
    pivot_longer(starts_with("All") | contains("FMT")) |>
    separate_wider_delim(
        name,
        delim = "_", names = c("cor_group", "param")
    ) |>
    pivot_wider(names_from = param, values_from = value)

# Only consider interactions where the entire model fits well enough.
library(ggplot2)
anansiLong <- anansiLong[anansiLong$full_p.values < 0.05, ]



ggplot(
    data = anansiLong,
    aes(
        x = r.values,
        y = feature_X,
        fill = cor_group,
        alpha = disjointed_Legend_p.values < 0.05
    )
) +

    # Make a vertical dashed red line at x = 0
    geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +

    # Points show  raw correlation coefficients
    geom_point(shape = 21, size = 3) +

    # facet per compound
    ggforce::facet_col(~feature_Y, space = "free", scales = "free_y") +
    scale_fill_manual(values = c(
        "Young yFMT" = "#2166ac",
        "Aged oFMT" = "#b2182b",
        "Aged yFMT" = "#ef8a62",
        "All" = "gray"
    )) +
    theme_bw()
#> Warning: Using alpha for a discrete variable is not advised.


# Using miaViz style function:

p <- plotAnansi(anansi_out,
    association.type = "disjointed",
    model.var = "Legend",
    fill_by = "group",
    signif.threshold = 0.05,
    x_lab = "Pearson's rho"
)
p <- p +
    scale_fill_manual(values = c(
        "Young yFMT" = "#2166ac",
        "Aged oFMT" = "#b2182b",
        "Aged yFMT" = "#ef8a62",
        "All" = "gray"
    )) +
    theme_bw()

p

```
