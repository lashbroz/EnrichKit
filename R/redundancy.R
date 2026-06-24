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
#'   `redundancy_reason` columns. All original pathway rows are retained so the
#'   output can be cross-referenced across analyses run on different platforms
#'   or interrogated gene sets.
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

#' Build a Pathway Consolidation Crosswalk
#'
#' Converts a redundancy- or cascade-annotated pathway result table into a
#' representative-to-all-pathways map. This is useful when the same pathway
#' database is analyzed across different platforms: the retained representative
#' can be used for display, while all original pathway names remain available
#' for cross-referencing across protein, RNA, phosphosite, or other analyses.
#'
#' @param results A pathway result table returned by
#'   [reduce_redundant_pathways()], [consolidate_pathways_by_gene_separation()],
#'   or [cascade_threshold_consolidation()]$results.
#' @param pathway_col Column containing original pathway names.
#' @param representative_col Column containing representative pathway names. If
#'   absent and `cascade_retained` is present, pathways are assigned to the most
#'   recent retained pathway in cascade order.
#' @param collapse String used to collapse pathway names.
#'
#' @return A data frame with one row per representative pathway and columns for
#'   all original pathway names assigned to that representative.
#' @export
pathway_consolidation_map <- function(results,
                                      pathway_col = "pathway",
                                      representative_col = "representative_pathway",
                                      collapse = ";") {
  if (!is.data.frame(results)) {
    stop("`results` must be a data frame.")
  }
  if (!pathway_col %in% colnames(results)) {
    stop("Missing pathway column: ", pathway_col)
  }

  pathways <- as.character(results[[pathway_col]])
  if (representative_col %in% colnames(results)) {
    representatives <- as.character(results[[representative_col]])
    representatives[is.na(representatives) | !nzchar(representatives)] <-
      pathways[is.na(representatives) | !nzchar(representatives)]
  } else if ("cascade_retained" %in% colnames(results)) {
    ord <- if ("cascade_order" %in% colnames(results)) order(results$cascade_order) else seq_len(nrow(results))
    representatives <- rep(NA_character_, length(pathways))
    current_representative <- NA_character_
    for (i in ord) {
      if (isTRUE(results$cascade_retained[i]) || is.na(current_representative)) {
        current_representative <- pathways[i]
      }
      representatives[i] <- current_representative
    }
  } else {
    stop(
      "Could not infer pathway representatives. Provide `",
      representative_col,
      "` or a cascade result with `cascade_retained`."
    )
  }

  groups <- split(pathways, representatives)
  rep_order <- unique(representatives)
  groups <- groups[rep_order]

  out <- data.frame(
    representative_pathway = names(groups),
    n_consolidated_pathways = lengths(groups),
    consolidated_pathways = vapply(groups, function(x) paste(unique(x), collapse = collapse), character(1)),
    stringsAsFactors = FALSE
  )

  out$dropped_pathways <- vapply(seq_along(groups), function(i) {
    dropped <- setdiff(unique(groups[[i]]), names(groups)[i])
    paste(dropped, collapse = collapse)
  }, character(1))

  if ("database" %in% colnames(results)) {
    idx <- match(out$representative_pathway, pathways)
    out$representative_database <- as.character(results$database[idx])
  }

  out
}

#' Cross-Reference Pathway Results Across Analyses
#'
#' Builds a pathway-key table across multiple enrichment analyses. Each original
#' pathway name is treated as the stable key, while optional representative
#' pathway columns preserve how that pathway was displayed after redundancy
#' reduction or cascade consolidation. This supports cross-platform comparisons,
#' for example asking whether a pathway seen in a protein analysis is also
#' present in an RNA, phosphosite, or collaborator-subset analysis.
#'
#' @param results_list Named list of pathway result tables.
#' @param pathway_col Column containing original pathway names.
#' @param representative_col Optional column containing consolidated/display
#'   pathway names. Missing columns are filled with the original pathway name.
#' @param p_col Optional p-value column.
#' @param fdr_col Optional FDR column.
#' @param direction_col Optional direction/effect column.
#' @param database_col Optional database/source label column.
#' @param keep_all_rows If `FALSE`, keep the strongest row per
#'   analysis/pathway by FDR then p-value.
#'
#' @return A list with `long` pathway-analysis records and a `wide` presence
#'   table keyed by original pathway name.
#' @export
cross_reference_pathways <- function(results_list,
                                     pathway_col = "pathway",
                                     representative_col = "representative_pathway",
                                     p_col = "p",
                                     fdr_col = "fdr",
                                     direction_col = "dir",
                                     database_col = "database",
                                     keep_all_rows = FALSE) {
  if (!is.list(results_list) || is.null(names(results_list))) {
    stop("`results_list` must be a named list of result tables.")
  }
  if (any(!nzchar(names(results_list)))) {
    stop("All entries in `results_list` must have non-empty names.")
  }

  long <- do.call(rbind, lapply(names(results_list), function(analysis) {
    x <- results_list[[analysis]]
    if (!is.data.frame(x)) {
      stop("Result entry `", analysis, "` is not a data frame.")
    }
    if (!pathway_col %in% colnames(x)) {
      stop("Result entry `", analysis, "` is missing pathway column: ", pathway_col)
    }

    pathway <- as.character(x[[pathway_col]])
    representative <- if (representative_col %in% colnames(x)) {
      as.character(x[[representative_col]])
    } else {
      pathway
    }
    representative[is.na(representative) | !nzchar(representative)] <-
      pathway[is.na(representative) | !nzchar(representative)]

    data.frame(
      analysis = analysis,
      pathway_key = pathway,
      representative_pathway = representative,
      database = get_optional_column(x, database_col, NA_character_),
      p = suppressWarnings(as.numeric(get_optional_column(x, p_col, NA_real_))),
      fdr = suppressWarnings(as.numeric(get_optional_column(x, fdr_col, NA_real_))),
      direction = suppressWarnings(as.numeric(get_optional_column(x, direction_col, NA_real_))),
      stringsAsFactors = FALSE
    )
  }))

  if (!isTRUE(keep_all_rows) && nrow(long) > 0) {
    ord <- order(
      long$analysis,
      long$pathway_key,
      ifelse(is.na(long$fdr), Inf, long$fdr),
      ifelse(is.na(long$p), Inf, long$p)
    )
    long <- long[ord, , drop = FALSE]
    long <- long[!duplicated(paste(long$analysis, long$pathway_key, sep = "\r")), , drop = FALSE]
  }
  rownames(long) <- NULL

  analyses <- names(results_list)
  keys <- sort(unique(long$pathway_key))
  wide <- data.frame(pathway_key = keys, stringsAsFactors = FALSE)
  for (analysis in analyses) {
    idx <- match(keys, long$pathway_key[long$analysis == analysis])
    sub <- long[long$analysis == analysis, , drop = FALSE]
    present <- !is.na(idx)
    wide[[paste0(analysis, "_present")]] <- present
    wide[[paste0(analysis, "_representative")]] <- ifelse(present, sub$representative_pathway[idx], NA_character_)
    wide[[paste0(analysis, "_fdr")]] <- ifelse(present, sub$fdr[idx], NA_real_)
    wide[[paste0(analysis, "_direction")]] <- ifelse(present, sub$direction[idx], NA_real_)
  }
  wide$n_analyses_present <- rowSums(wide[paste0(analyses, "_present")])
  wide <- wide[order(-wide$n_analyses_present, wide$pathway_key), , drop = FALSE]
  rownames(wide) <- NULL

  list(long = long, wide = wide)
}

get_optional_column <- function(x, col, default) {
  if (col %in% colnames(x)) {
    return(x[[col]])
  }
  rep(default, nrow(x))
}

two_set_similarity <- function(a, b) {
  intersection <- intersect(a, b)
  union_set <- union(a, b)
  list(
    jaccard = if (length(union_set) == 0) NA_real_ else length(intersection) / length(union_set),
    gene_separation = length(setdiff(union_set, intersection))
  )
}
