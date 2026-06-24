#' Pairwise Pathway Similarity
#'
#' Calculates pairwise Jaccard similarity and gene separation between pathway
#' gene sets.
#'
#' @param gene_sets Named list of gene sets or an `EnrichKit_pathway_db` object.
#'
#' @return A data frame of pairwise pathway similarities.
#' @export
pathway_similarity <- function(gene_sets) {
  sets <- as_gene_sets(gene_sets)
  if (length(sets) < 2) {
    return(data.frame())
  }
  pairs <- utils::combn(names(sets), 2, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(pair) {
    a <- sets[[pair[1]]]
    b <- sets[[pair[2]]]
    intersection <- intersect(a, b)
    union <- union(a, b)
    data.frame(
      pathway_a = pair[1],
      pathway_b = pair[2],
      n_a = length(a),
      n_b = length(b),
      n_intersection = length(intersection),
      n_union = length(union),
      jaccard = if (length(union) == 0) NA_real_ else length(intersection) / length(union),
      gene_separation = length(setdiff(union, intersection)),
      stringsAsFactors = FALSE
    )
  }))
}

#' Reduce Redundant Pathways
#'
#' Greedy redundancy reduction. Pathways are visited in result-table order, and
#' later pathways are dropped when they overlap a retained pathway above the
#' chosen Jaccard cutoff or at/below the chosen gene-separation cutoff.
#'
#' @param results Enrichment result table.
#' @param gene_sets Named list of gene sets or an `EnrichKit_pathway_db` object.
#' @param pathway_col Result column containing pathway names.
#' @param jaccard_cutoff Drop pathways with Jaccard similarity greater than or
#'   equal to this cutoff. Use `NULL` to disable.
#' @param gene_separation_cutoff Drop pathways with gene separation less than or
#'   equal to this cutoff. Use `NULL` to disable.
#'
#' @return `results` with `redundant`, `representative_pathway`, and
#'   `redundancy_reason` columns.
#' @export
reduce_redundant_pathways <- function(results,
                                      gene_sets,
                                      pathway_col = "pathway",
                                      jaccard_cutoff = 0.8,
                                      gene_separation_cutoff = NULL) {
  if (!pathway_col %in% colnames(results)) {
    stop("Missing pathway column: ", pathway_col)
  }
  sets <- as_gene_sets(gene_sets)
  pathways <- as.character(results[[pathway_col]])
  missing_sets <- setdiff(pathways, names(sets))
  if (length(missing_sets) > 0) {
    warning(length(missing_sets), " result pathway(s) missing from gene sets.")
  }

  retained <- character(0)
  redundant <- rep(FALSE, length(pathways))
  representative <- rep(NA_character_, length(pathways))
  reason <- rep(NA_character_, length(pathways))

  for (i in seq_along(pathways)) {
    p <- pathways[i]
    if (!p %in% names(sets)) {
      next
    }
    matched <- FALSE
    for (r in retained) {
      sim <- two_set_similarity(sets[[p]], sets[[r]])
      drop_jaccard <- !is.null(jaccard_cutoff) && is.finite(sim$jaccard) && sim$jaccard >= jaccard_cutoff
      drop_sep <- !is.null(gene_separation_cutoff) && sim$gene_separation <= gene_separation_cutoff
      if (drop_jaccard || drop_sep) {
        redundant[i] <- TRUE
        representative[i] <- r
        reason[i] <- if (drop_jaccard) "jaccard" else "gene_separation"
        matched <- TRUE
        break
      }
    }
    if (!matched) {
      retained <- c(retained, p)
      representative[i] <- p
      reason[i] <- "retained"
    }
  }

  results$redundant <- redundant
  results$representative_pathway <- representative
  results$redundancy_reason <- reason
  results
}

#' Consolidate Pathways by Gene-Separation Threshold
#'
#' Convenience wrapper for fixed gene-separation redundancy reduction, matching
#' the manuscript workflow where closely related pathways are consolidated by a
#' stable separation threshold.
#'
#' @inheritParams reduce_redundant_pathways
#' @param gene_separation_threshold Maximum separation for consolidation.
#'
#' @return A redundancy-annotated result table.
#' @export
consolidate_pathways_by_gene_separation <- function(results,
                                                    gene_sets,
                                                    pathway_col = "pathway",
                                                    gene_separation_threshold = 10) {
  reduce_redundant_pathways(
    results = results,
    gene_sets = gene_sets,
    pathway_col = pathway_col,
    jaccard_cutoff = NULL,
    gene_separation_cutoff = gene_separation_threshold
  )
}

two_set_similarity <- function(a, b) {
  intersection <- intersect(a, b)
  union_set <- union(a, b)
  list(
    jaccard = if (length(union_set) == 0) NA_real_ else length(intersection) / length(union_set),
    gene_separation = length(setdiff(union_set, intersection))
  )
}
