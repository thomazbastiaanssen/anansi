# Methods for MultiFactor S7 container class

- [`droplevels()`](https://rdrr.io/r/base/droplevels.html) :

  `droplevels(MultiFactor)` returns a `MultiFactor` with unused levels
  removed, Analogous to the `factor` method.

## Arguments

- ...:

  `i,j` indices specifying elements to extract or replace. Indices are
  numeric or character vectors or empty (missing) or NULL. Numeric
  values are coerced to integer or whole numbers as by as.integer or for
  large values by trunc (and hence truncated towards zero). Character
  vectors will be matched to the names of the object.

- value:

  Replacement value, typically of same type as that which is to be
  replaced.

- exclude:

  `NULL` or `Named character list` of similar structure as
  `levels(MultiFactor)`. Which levels to drop from output.

- select:

  `NULL` or `Named character list` of similar structure as
  `levels(MultiFactor)`. Which levels to keep in output.

- use.names, ignore.mcols:

  For compatibility, not used.

- x:

  `MultiFactor`

- format:

  `Character scalar`, controls output format by package name. `"igraph"`
  and `"graph"` are supported.

## Value

A MultiFactor

a specified graph object.

## Details

Only one of `select` and `exclude` should be provided, as they are each
others complement.

## See also

[`igraph::graph_from_data_frame()`](https://r.igraph.org/reference/graph_from_data_frame.html)
and
[`igraph::as_graphnel()`](https://r.igraph.org/reference/as_graphnel.html),
which are used under the hood, from
[`igraph::igraph()`](https://r.igraph.org/reference/aaa-igraph-package.html)
package.

## Examples

``` r
# Setup
x <- MultiFactor(kegg_link())
x
#> An anansi::MultiFactor S7_object,
#>     3 feature types across 2 edge lists.
#> 
#>          ec   ko  cpd
#> ec2ko  5581 8817    .
#> ec2cpd 5934    . 7802
#> 
#> Values represent unique feature names in that edge list.
#> 
#> Levels:
#> 
#> ec  : 6690 Levels: 1.1.1.1 1.1.1.2 ... 5.3.1.37 
#> ko  : 8817 Levels: K00001 K00002 ... K28089 
#> cpd : 7802 Levels: C00001 C00002 ... C23000 

# Basic properties
dim(x)
#> [1] 2 3
dimnames(x)
#> [[1]]
#> [1] "ec2ko"  "ec2cpd"
#> 
#> [[2]]
#> [1] "ec"  "ko"  "cpd"
#> 

# Factor-like properties
head(levels(x)$ko)
#> [1] "K00001" "K00002" "K00003" "K00004" "K00005" "K00006"
droplevels(x)
#> An anansi::MultiFactor S7_object,
#>     3 feature types across 2 edge lists.
#> 
#>          ec   ko  cpd
#> ec2ko  5581 8817    .
#> ec2cpd 5934    . 7802
#> 
#> Values represent unique feature names in that edge list.
#> 
#> Levels:
#> 
#> ec  : 6690 Levels: 1.1.1.1 1.1.1.2 ... 5.3.1.37 
#> ko  : 8817 Levels: K00001 K00002 ... K28089 
#> cpd : 7802 Levels: C00001 C00002 ... C23000 
head(unfactor(x)$ec2ko)
#>          ec     ko
#> 1   1.1.1.1 K00001
#> 2   1.1.1.2 K00002
#> 3   1.1.1.3 K00003
#> 4   1.1.1.4 K00004
#> 5 1.1.1.303 K00004
#> 6   1.1.1.6 K00005

# Extract common output formats
getEdgeList(x)
#>        V1  V2
#> ec2ko  ec  ko
#> ec2cpd ec cpd

droplevels(x, exclude = list(ko = "K00001"))
#> An anansi::MultiFactor S7_object,
#>     3 feature types across 2 edge lists.
#> 
#>          ec   ko  cpd
#> ec2ko  5581 8816    .
#> ec2cpd 5934    . 7802
#> 
#> Values represent unique feature names in that edge list.
#> 
#> Levels:
#> 
#> ec  : 6690 Levels: 1.1.1.1 1.1.1.2 ... 5.3.1.37 
#> ko  : 8816 Levels: K00002 K00003 ... K28089 
#> cpd : 7802 Levels: C00001 C00002 ... C23000 
droplevels(x, select = list(ko = "K00001"))
#> An anansi::MultiFactor S7_object,
#>     3 feature types across 2 edge lists.
#> 
#>          ec ko  cpd
#> ec2ko     1  1    .
#> ec2cpd 5934  . 7802
#> 
#> Values represent unique feature names in that edge list.
#> 
#> Levels:
#> 
#> ec  : 5934 Levels: 1.1.1.1 1.1.1.2 ... 5.3.1.37 
#> ko  : 1 Levels: K00001 
#> cpd : 7802 Levels: C00001 C00002 ... C23000 
# Generate an igraph object
g <- getGraph(x = kegg_link(), format = "igraph")
plot(g)

```
