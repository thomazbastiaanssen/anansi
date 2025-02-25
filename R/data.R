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

#' @rdname kegg_link
#' @format \code{ec2ko}: a \code{data.frame} of two columns, named \code{"ec"} 
#' and \code{"ko"}. 
#' The IDs refer to KEGG orthologues. Enzyme commission numbers, ecs, typically
#' describe reactions captured by them.
#'
#' @source \code{ec2ko}: Adapted from \url{https://www.genome.jp/kegg/}, using 
#' \code{KEGGREST}. Script to generate available in example. 
#'                     
"ec2ko"

#' @rdname kegg_link
#' @format \code{ec2cpd}: a \code{data.frame} of two columns, named \code{"ec"} 
#' and \code{"cpd"}. 
#' The IDs refer to compounds in the KEGG database. Enzyme commission numbers, 
#' ecs, typically describe reactions either producing or requiring them.
#'
#' @source \code{ec2cpd}: Adapted from \url{https://www.genome.jp/kegg/} using 
#' \code{KEGGREST}. Script to generate available in example. 
#'
"ec2cpd"

#' Use linking data from the KEGG database. 
#' @description 
#' \code{kegg_links} is a convenience function to return a list 
#' containing two \code{data.frame}s; \code{ec2cpd} and \code{ec2ko}. This will 
#' be their most likely use.
#'  
#' \code{ec2cpd} and \code{ec2ko} are two \code{data.frame}s, used 
#' to link ko, ecs and cpd identifiers in the KEGG database.  
#' 

#' @returns \code{kegg_link} returns a list containing the two aforementioned 
#' data.frames, \code{ec2cpd} and \code{ec2ko}.
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
kegg_link <- function() list(anansi::ec2ko, anansi::ec2cpd)

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
