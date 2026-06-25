#' Get the Packaged Kids First gosets.all Database
#'
#' Loads and optionally filters the packaged Kids First `gosets.all` object.
#' This is the easiest way to use the current Kids First pathway database shipped
#' with EnrichKit.
#'
#' @param gene Optional gene universe. If supplied, each pathway is intersected
#'   with these genes before size filtering.
#' @param min_genes Minimum retained pathway size.
#' @param max_genes Maximum retained pathway size.
#' @param deduplicate If `TRUE`, remove duplicate gene sets after filtering.
#' @param as_pathway_db If `TRUE`, return an `EnrichKit_pathway_db` object with
#'   source labels from `kfirst_gosets_source`. If `FALSE`, return a named list.
#'
#' @return A named list of gene sets or an `EnrichKit_pathway_db` object.
#' @export
get_kfirst_gosets <- function(gene = NULL,
                              min_genes = 5,
                              max_genes = 250,
                              deduplicate = FALSE,
                              as_pathway_db = FALSE) {
  data_env <- new.env(parent = emptyenv())
  utils::data("kfirst_gosets_all", package = "EnrichKit", envir = data_env)
  utils::data("kfirst_gosets_source", package = "EnrichKit", envir = data_env)
  kfirst_gosets_all <- get("kfirst_gosets_all", envir = data_env)
  kfirst_gosets_source <- get("kfirst_gosets_source", envir = data_env)

  out <- kfirst_gosets_all
  if (!is.null(gene)) {
    gene <- unique(as.character(gene))
    out <- lapply(out, function(x) sort(unique(intersect(x, gene))))
  }

  n_genes <- lengths(out)
  out <- out[n_genes >= min_genes & n_genes <= max_genes]

  if (deduplicate) {
    collapsed <- vapply(out, paste, character(1), collapse = ",")
    out <- out[match(unique(collapsed), collapsed)]
  }

  if (isTRUE(as_pathway_db)) {
    idx <- match(names(out), kfirst_gosets_source$pathway)
    database <- kfirst_gosets_source$source[idx]
    database[is.na(database)] <- "KidsFirst_gosets_all"
    out <- make_pathway_db(
      out,
      database = stats::setNames(database, names(out)),
      source = "KidsFirst_packaged_gosets_all",
      min_size = min_genes,
      max_size = max_genes,
      deduplicate_identical = FALSE
    )
  }

  out
}
