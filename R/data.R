#' Current Kids First Pathway Database
#'
#' The current Kids First pathway database used by EnrichKit.
#'
#' The packaged database was filtered to the Kids First measured assay
#' background using pathway-size limits of 5 to 250 genes.
#'
#' @format A named list. Each element is a character vector of gene symbols for
#'   one pathway/gene set.
#' @source Packaged Kids First pathway database.
"kfirst_gosets_all"

#' Kids First Pathway Source Table
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

#' Kids First Measured Assay Background
#'
#' Measured assay background used to filter [kfirst_gosets_all].
#'
#' @format Character vector of gene symbols.
"kfirst_gene_universe"

#' Kids First Pathway Metadata
#'
#' Metadata describing the packaged [kfirst_gosets_all] object.
#'
#' @format A list with source labels, pathway-size filters, pathway counts, and
#'   the size of the Kids First measured assay background.
"kfirst_gosets_metadata"
