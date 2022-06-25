#' Dictionary to link KEGG compounds to KEGG orthologues
#'
#' @format a list of 18844 named vectors. The names refer to compounds in the KEGG database.
#' The entries in the vectors refer to the KEGG orthologues that are associated to those compounds;
#' typically either producing or requiring them.
#'
#' @source Adapted from \url{https://www.genome.jp/kegg/}
#'
"anansi_dic"

#' Dictionary to link KEGG orthologues based on their pathway membership.
#'
#' @format a list of 467 named vectors. The names refer to pathways in the KEGG database.
#' The entries in the vectors refer to the KEGG orthologues that are make up those pathways.
#'
#' @source Adapted from \url{https://www.genome.jp/kegg/}
#'
"KO_member_dic"

#' Lookup table to convert KEGG compound IDs to their human-readable counterparts
#'
#' @format A data.frame with two columns and 18871 rows. The first column refers to all compounds in the KEGG database.
#' The second column refers to their human readable counterpart.
#' @source Adapted from \url{https://www.genome.jp/kegg/}
#'
"cpd_translation"

#' Lookup table to convert KEGG orthologue (KO) IDs to their human-readable counterparts
#'
#' @format A data.frame with two columns and 24620 rows. The first column refers to all KEGG orthologues (KOs) in the KEGG database.
#' The second column refers to their human readable counterpart.
#' @source Adapted from \url{https://www.genome.jp/kegg/}
#'
"KO_translation"

#' Snippet of the CLR-transformed hippocampal metabolomics data from the FMT Aging study.
#'
#' @format A matrix object with three rows, compounds, and 36 columns, samples.
#' @source \url{https://doi.org/10.1038/s43587-021-00093-9}
#'
"FMT_metab"

#' Snippet of the untransformed inferred functional data from the FMT Aging study.
#'
#' @description Piphillin was used to infer funcitons from the 16S sequencing data in terms of KOs.
#' Unfortunately, the Piphillin algorithm is proprietary and has since been taken down.
#'
#' @format A marix object with 6474 rows, KOs, and 36 columns, samples.
#' @source \url{https://doi.org/10.1038/s43587-021-00093-9}
#'
"FMT_KOs"

#' Snippet of the metadata from the FMT Aging study.
#'
#' @description There were three treatment groups in the study that all received faecal microbiota transplantation (FMT).
#' Young mice that received FMT from young mice (Young yFMT),
#' aged mice that received FMT from aged mice (Aged oFMT) and aged mice that received FMT from young mice (Aged yFMT).
#'
#' @format A data.frame object with 36 rows, samples, and two columns, denoting sample ID and treatemtn group, respectively.
#' @source \url{https://doi.org/10.1038/s43587-021-00093-9}
#'
"FMT_metadata"
