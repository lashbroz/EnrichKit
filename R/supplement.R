#' Format Publication Supplementary Enrichment Table
#'
#' Produces a stable, manuscript-friendly enrichment table with consistent
#' column names, database labels, p-values, adjusted p-values, effect direction,
#' gene counts, background size, and optional redundancy/module metadata.
#'
#' @param results Enrichment result data frame.
#' @param analysis_name Analysis label to add to the table.
#' @param pathway_col Pathway column.
#' @param database_col Database column.
#' @param p_col P-value column.
#' @param fdr_col FDR column.
#' @param direction_col Direction/effect sign column.
#' @param effect_col Optional effect-size column.
#' @param sort Sort by FDR then p-value.
#'
#' @return A data frame with stable supplementary-table columns.
#' @export
format_enrichment_supplement <- function(results,
                                         analysis_name = NA_character_,
                                         pathway_col = "pathway",
                                         database_col = "database",
                                         p_col = "p",
                                         fdr_col = "fdr",
                                         direction_col = "dir",
                                         effect_col = NULL,
                                         sort = TRUE) {
  required <- c(pathway_col, p_col, fdr_col)
  missing <- setdiff(required, colnames(results))
  if (length(missing) > 0) {
    stop("Missing required result columns: ", paste(missing, collapse = ", "))
  }
  if (!database_col %in% colnames(results)) {
    warning("Database column missing; using `unknown`.")
    results[[database_col]] <- "unknown"
  }
  if (!direction_col %in% colnames(results)) {
    results[[direction_col]] <- NA_real_
  }
  if (is.null(effect_col)) {
    effect_col <- if ("effect" %in% colnames(results)) "effect" else if ("odds_ratio" %in% colnames(results)) "odds_ratio" else NA_character_
  }

  out <- data.frame(
    analysis = analysis_name,
    pathway = as.character(results[[pathway_col]]),
    database = as.character(results[[database_col]]),
    method = if ("method" %in% colnames(results)) as.character(results$method) else NA_character_,
    alternative = if ("alternative" %in% colnames(results)) as.character(results$alternative) else NA_character_,
    direction = as.numeric(results[[direction_col]]),
    effect = if (!is.na(effect_col) && effect_col %in% colnames(results)) as.numeric(results[[effect_col]]) else NA_real_,
    p_value = as.numeric(results[[p_col]]),
    fdr = as.numeric(results[[fdr_col]]),
    fdr_scope = if ("fdr_scope" %in% colnames(results)) as.character(results$fdr_scope) else NA_character_,
    n_background = if ("n_background" %in% colnames(results)) as.integer(results$n_background) else NA_integer_,
    n_pathway = if ("n_pathway" %in% colnames(results)) as.integer(results$n_pathway) else NA_integer_,
    n_overlap = if ("n_overlap" %in% colnames(results)) as.integer(results$n_overlap) else NA_integer_,
    genes = if ("overlap_genes" %in% colnames(results)) as.character(results$overlap_genes) else if ("pathway_genes" %in% colnames(results)) as.character(results$pathway_genes) else NA_character_,
    stringsAsFactors = FALSE
  )

  optional <- intersect(
    c("redundant", "representative_pathway", "redundancy_reason", "module"),
    colnames(results)
  )
  for (col in optional) {
    out[[col]] <- results[[col]]
  }

  if (isTRUE(sort)) {
    out <- out[order(out$fdr, out$p_value), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' Write Publication Supplementary Enrichment Table
#'
#' @inheritParams format_enrichment_supplement
#' @param file Output TSV file path.
#' @param ... Additional arguments passed to [format_enrichment_supplement()].
#'
#' @return Invisibly returns `file`.
#' @export
write_enrichment_supplement <- function(results,
                                        file,
                                        analysis_name = NA_character_,
                                        ...) {
  out <- format_enrichment_supplement(results, analysis_name = analysis_name, ...)
  utils::write.table(out, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(file)
}
