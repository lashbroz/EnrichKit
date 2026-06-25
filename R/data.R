#' Current Kids First gosets.all Pathway Database
#'
#' The current Kids First `gosets.all` pathway database used by EnrichKit.
#'
#' This database was built from the HOPE pathway database after removing
#' `KEGG_MEDICUS*` pathways, supplemented with canonical KEGG pathways from
#' MSigDB C2 canonical pathways that were not already present. The combined
#' database was then filtered to the Kids First interrogated gene universe using
#' pathway-size limits of 5 to 250 genes.
#'
#' @format A named list. Each element is a character vector of gene symbols for
#'   one pathway/gene set.
#' @source HOPE pathway database plus canonical MSigDB C2 KEGG pathways.
"kfirst_gosets_all"

#' Kids First gosets.all Source Table
#'
#' Source/provenance table for [kfirst_gosets_all].
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{pathway}{Pathway/gene-set name.}
#'   \item{source}{Source database label.}
#'   \item{n_genes}{Number of retained genes after Kids First filtering.}
#' }
"kfirst_gosets_source"

#' Kids First Interrogated Gene Universe
#'
#' Gene universe used to filter [kfirst_gosets_all].
#'
#' @format Character vector of gene symbols.
"kfirst_gene_universe"

#' Kids First gosets.all Metadata
#'
#' Metadata describing the packaged [kfirst_gosets_all] object.
#'
#' @format A list with source labels, pathway-size filters, pathway counts, and
#'   the size of the Kids First gene universe.
"kfirst_gosets_metadata"
