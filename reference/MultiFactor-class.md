# MultiFactor S7 container class

`MultiFactor` is an S7 class to organize and manage multiple sets of
factors, for instance when tracing or converting feature IDs across
databases. Methods for `MultiFactor` aim to follow `factor` behaviour.

## Usage

``` r
MultiFactor(x, levels = NULL, drop.unmatched = FALSE)

## Constructor for `MultiFactor` objects
MultiFactor(x, levels = NULL, drop.unmatched = FALSE)
```

## Arguments

- x:

  `MultiFactor`

- levels:

  an optional named list of vectors of the unique values (as character
  strings) that x might have taken. The default is the unique set of
  values taken by lapply(x, as.character), sorted into increasing order
  of x.

- drop.unmatched:

  `Logical scalar` If `TRUE` (Default), for feature types that are seen
  at least twice, exclude features that only present in one of their
  respective link data frames.

## Value

a MultiFactor object.

## Details

The most straightforward way to construct a `MultiFactor` object is as a
named list of named data.frames. The columns of the data.frames indicate
the category of factor in that column.

A `MultiFactor` object presents itself similar to a `data.frame`, in the
sense that level types can be called as columns and individual
data.frame components can be called as rows.

## Slots

- `index`:

  Named `list` of named integer data frames of at least two columns
  each. The column names correspond to names in the `levels` slot.
  Similar to `factor`s, the integers in those columns correspond to the
  characters in that level. Accessed through regular list methods (e.g.,
  `[`, `[[`).

- `levels`:

  `Named list of character vectors`. Accessed through `levels(x)`

- `map`:

  `(sparse) Matrix` specifying which elements contain which levels.
  Accesses through `x@dictionary`.

## See also

[`MultiFactor-methods()`](https://thomazbastiaanssen.github.io/anansi/reference/MultiFactor-methods.md)

- [`kegg_link()`](https://thomazbastiaanssen.github.io/anansi/reference/kegg_link.md):
  for an example of valid input.

## Examples

``` r
# Generate some random linkage input
x <- data.frame(
    a = sample(letters[seq(3)], 10, replace = TRUE),
    A = sample(LETTERS[seq(3)], 10, replace = TRUE)
)
MultiFactor(x)
#> An anansi::MultiFactor S7_object,
#>     2 feature types across 1 edge lists.
#> 
#>   a A
#> x 3 2
#> 
#> Values represent unique feature names in that edge list.
#> 
#> Levels:
#> 
#> a : 3 Levels: a b c 
#> A : 2 Levels: C B 
```
