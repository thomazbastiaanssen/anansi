#' @rdname kegg_link
#' @format `ec2ko`: a `data.frame` of two columns, named `"ec"`
#' and `"ko"`.
#' The IDs refer to KEGG orthologues. Enzyme commission numbers, ecs, typically
#' describe reactions captured by them.
#'
#' @source `ec2ko`: Adapted from <https://www.genome.jp/kegg/>, using
#' `KEGGREST`. Script to generate available in example.
#'
"ec2ko"

#' @rdname kegg_link
#' @format `ec2cpd`: a `data.frame` of two columns, named `"ec"`
#' and `"cpd"`.
#' The IDs refer to compounds in the KEGG database. Enzyme commission numbers,
#' ecs, typically describe reactions either producing or requiring them.
#'
#' @source `ec2cpd`: Adapted from <https://www.genome.jp/kegg/> using
#' `KEGGREST`. Script to generate available in example.
#'
"ec2cpd"

#' Use linking data from the KEGG database.
#' @description
#' `kegg_links` is a convenience function to return a list
#' containing two `data.frame`s; `ec2cpd` and `ec2ko`. This will
#' be their most likely use.
#'
#' `ec2cpd` and `ec2ko` are two `data.frame`s, used
#' to link ko, ecs and cpd identifiers in the KEGG database.
#'

#' @returns `kegg_link` returns a list containing the two aforementioned
#' data.frames, `ec2cpd` and `ec2ko`.
#' @examples
#' kegg_link()
#'
#' # Generate ec2ko and ec2cpd:
#' # Don't download during tests. set to `TRUE` to download.
#' dry_run <- TRUE
#'
#' if(!dry_run) {
#'
#' ec2ko <- KEGGREST::keggLink("ec", "ko")
#' ec2ko <- data.frame(ec = gsub("ec:","",  x = ec2ko, fixed = TRUE),
#'                     ko = gsub("ko:", "", x = names(ec2ko), fixed = TRUE),
#'                     row.names = NULL)
#'
#' ec2cpd <- KEGGREST::keggLink("ec", "cpd")
#' ec2cpd <- data.frame(ec  = gsub("ec:","",   x = ec2cpd, fixed = TRUE),
#'                      cpd = gsub("cpd:", "", x = names(ec2cpd), fixed = TRUE),
#'                      row.names = NULL)
#' }
#'
#' @export
#'
kegg_link <- function() list(ec2ko = anansi::ec2ko, ec2cpd = anansi::ec2cpd)

#' Snippet of the CLR-transformed hippocampal metabolomics data from the FMT Aging study.
#'
#' @format A matrix object with three rows, compounds, and 36 columns, samples.
#' @source \doi{10.1038/s43587-021-00093-9}
#'
"FMT_metab"

#' Snippet of the untransformed inferred functional data from the FMT Aging study.
#'
#' @description Piphillin was used to infer functions from the 16S sequencing data in terms of KOs.
#' Unfortunately, the Piphillin algorithm is proprietary and has since been taken down.
#'
#' @format A marix object with 6474 rows, KOs, and 36 columns, samples.
#' @source \doi{10.1038/s43587-021-00093-9}
#'
"FMT_KOs"

#' Snippet of the metadata from the FMT Aging study.
#'
#' @description There were three treatment groups in the study that all received faecal microbiota transplantation (FMT).
#' Young mice that received FMT from young mice (Young yFMT),
#' aged mice that received FMT from aged mice (Aged oFMT) and aged mice that received FMT from young mice (Aged yFMT).
#'
#' @format A data.frame object with 36 rows, samples, and two columns, denoting sample ID and treatemtn group, respectively.
#' @source \doi{10.1038/s43587-021-00093-9}
#'
"FMT_metadata"
