#' Clean a List of Feature Sets
#'
#' Removes duplicated features within each set, intersects each set with an
#' optional universe, and filters sets by size.
#'
#' @param feature_sets Named list of character vectors.
#' @param universe Optional character vector defining the eligible feature universe.
#' @param min_size Minimum retained set size.
#' @param max_size Maximum retained set size.
#'
#' @return A named list of cleaned feature sets.
#' @export
clean_feature_sets <- function(feature_sets,
                               universe = NULL,
                               min_size = 1,
                               max_size = Inf) {
  clean_gene_sets(
    gene_sets = feature_sets,
    universe = universe,
    min_size = min_size,
    max_size = max_size
  )
}

#' Rank Features by a Numeric Statistic
#'
#' Converts a named numeric vector or two-column data frame into a descending
#' feature ranking.
#'
#' @param x Named numeric vector or data frame.
#' @param feature_col Feature column name when `x` is a data frame.
#' @param value_col Numeric statistic column name when `x` is a data frame.
#' @param decreasing Sort larger values first.
#'
#' @return A named numeric vector ordered by rank.
#' @export
rank_features <- function(x,
                          feature_col = "feature",
                          value_col = "statistic",
                          decreasing = TRUE) {
  if (is.data.frame(x)) {
    if (!all(c(feature_col, value_col) %in% colnames(x))) {
      stop("Data frame must contain `feature_col` and `value_col`.")
    }
    values <- x[[value_col]]
    names(values) <- as.character(x[[feature_col]])
  } else {
    values <- x
  }

  if (is.null(names(values))) {
    stop("Ranking vector must be named.")
  }

  values <- as.numeric(values)
  values <- values[!is.na(names(values)) & nzchar(names(values)) & is.finite(values)]
  values[order(values, decreasing = decreasing)]
}

#' Over-Representation Analysis
#'
#' Performs one-sided Fisher enrichment for a hit list against named feature sets.
#'
#' @param hits Character vector of selected features.
#' @param feature_sets Named list of feature sets.
#' @param universe Character vector of eligible background features. If omitted,
#'   the union of hits and feature sets is used.
#' @param min_size Minimum set size after universe filtering.
#' @param max_size Maximum set size after universe filtering.
#' @param p_adjust_method P-value adjustment method passed to [stats::p.adjust()].
#' @param fdr_by_database If `TRUE`, adjust FDR separately within each database
#'   label. If `FALSE`, adjust across all returned pathways together.
#' @param alternative Fisher alternative. `"greater"` tests over-representation
#'   of hits in the pathway, `"less"` tests under-representation, and
#'   `"two.sided"` tests either direction.
#'
#' @return A data frame with enrichment statistics.
#' @export
calc_ora <- function(hits,
                     feature_sets,
                     universe = NULL,
                     min_size = 1,
                     max_size = Inf,
                     p_adjust_method = "BH",
                     fdr_by_database = TRUE,
                     alternative = c("greater", "less", "two.sided")) {
  fisher_enrichment(
    hits = hits,
    pathway_db = feature_sets,
    background = universe,
    min_size = min_size,
    max_size = max_size,
    p_adjust_method = p_adjust_method,
    fdr_by_database = fdr_by_database,
    alternative = alternative
  )
}

#' Fisher/Exact-Test Pathway Enrichment
#'
#' Performs one-sided Fisher enrichment for a hit list against pathway gene sets.
#' This is the package version of the exact-test enrichment pattern used in the
#' HOPE/CPTAC scripts. The background should be the analysis-specific measured
#' feature universe, not all genes in the genome.
#'
#' @param hits Character vector of selected features.
#' @param pathway_db An `EnrichKit_pathway_db` object or named list of gene sets.
#' @param background Character vector of eligible background features. If
#'   omitted, the union of hits and gene sets is used.
#' @param min_size Minimum pathway size after background matching.
#' @param max_size Maximum pathway size after background matching.
#' @param p_adjust_method P-value adjustment method.
#' @param fdr_by_database If `TRUE`, adjust FDR separately within each database
#'   label. If `FALSE`, adjust across all returned pathways together.
#' @param alternative Fisher alternative. `"greater"` tests pathway
#'   over-representation among hits; `"less"` tests under-representation;
#'   `"two.sided"` performs a two-sided exact test.
#'
#' @return A data frame with pathway-level Fisher enrichment statistics.
#' @export
fisher_enrichment <- function(hits,
                              pathway_db,
                              background = NULL,
                              min_size = 1,
                              max_size = Inf,
                              p_adjust_method = "BH",
                              fdr_by_database = TRUE,
                              alternative = c("greater", "less", "two.sided")) {
  alternative <- match.arg(alternative)
  hits <- unique(as.character(hits))
  hits <- hits[!is.na(hits) & nzchar(hits)]
  gene_sets <- as_gene_sets(pathway_db)

  if (is.null(background)) {
    background <- unique(c(hits, unlist(gene_sets, use.names = FALSE)))
  }
  background <- unique(as.character(background))
  background <- background[!is.na(background) & nzchar(background)]

  hits <- intersect(hits, background)
  gene_sets <- clean_gene_sets(
    gene_sets = gene_sets,
    universe = background,
    min_size = min_size,
    max_size = max_size
  )

  if (length(gene_sets) == 0) {
    return(data.frame())
  }

  meta <- get_pathway_meta_for_results(pathway_db, names(gene_sets))

  res <- lapply(names(gene_sets), function(set_name) {
    set <- gene_sets[[set_name]]
    overlap <- intersect(hits, set)

    a <- length(overlap)
    b <- length(hits) - a
    c <- length(set) - a
    d <- length(background) - a - b - c

    contingency <- matrix(
      c(a, b, c, d),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        hit_status = c("hit", "not_hit"),
        pathway_status = c("in_pathway", "not_in_pathway")
      )
    )

    ft <- stats::fisher.test(contingency, alternative = alternative)

    data.frame(
      pathway = set_name,
      database = meta$database[match(set_name, meta$pathway)],
      n_background = length(background),
      n_hits = length(hits),
      n_pathway = length(set),
      n_overlap = a,
      odds_ratio = unname(ft$estimate),
      p = ft$p.value,
      dir = sign(unname(ft$estimate) - 1),
      alternative = alternative,
      overlap_genes = paste(overlap, collapse = ";"),
      method = "fisher",
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, res)
  if (isTRUE(fdr_by_database)) {
    out <- adjust_fdr_by_database(out, p_col = "p", database_col = "database", method = p_adjust_method)
    out$fdr_scope <- "database"
  } else {
    out$fdr <- stats::p.adjust(out$p, method = p_adjust_method)
    out$fdr_scope <- "global"
  }
  out$signed.p <- signed_log10_p(out$p, out$dir)
  out$signed.fdr <- signed_log10_p(out$fdr, out$dir)
  out[order(out$p, -out$n_overlap), ]
}

#' Wilcoxon Pathway Enrichment
#'
#' Tests whether pathway-member feature scores differ from background
#' non-member scores using a Wilcoxon rank-sum test. Direction is based on the
#' difference between pathway-member and non-member medians.
#'
#' @param feature_scores Named numeric vector of feature-level scores.
#' @param pathway_db An `EnrichKit_pathway_db` object or named list.
#' @param background Optional analysis-specific background features.
#' @param alternative Wilcoxon alternative.
#'   `"greater"` tests whether pathway-member scores are larger than
#'   non-member background scores; `"less"` tests whether they are smaller;
#'   `"two.sided"` tests either shift.
#' @param min_size Minimum pathway size after background matching.
#' @param max_size Maximum pathway size after background matching.
#' @param p_adjust_method P-value adjustment method.
#' @param fdr_by_database If `TRUE`, adjust FDR separately within each database
#'   label. If `FALSE`, adjust across all returned pathways together.
#'
#' @return A data frame of pathway-level Wilcoxon enrichment results.
#' @export
wilcox_enrichment <- function(feature_scores,
                              pathway_db,
                              background = NULL,
                              alternative = c("two.sided", "greater", "less"),
                              min_size = 1,
                              max_size = Inf,
                              p_adjust_method = "BH",
                              fdr_by_database = TRUE) {
  alternative <- match.arg(alternative)
  feature_scores <- clean_named_numeric(feature_scores)
  if (is.null(background)) {
    background <- names(feature_scores)
  }
  background <- intersect(unique(as.character(background)), names(feature_scores))
  if (length(background) < 2) {
    stop("Background must contain at least two scored features.")
  }

  sets <- clean_gene_sets(as_gene_sets(pathway_db), universe = background, min_size = min_size, max_size = max_size)
  if (length(sets) == 0) {
    return(data.frame())
  }
  meta <- get_pathway_meta_for_results(pathway_db, names(sets))

  res <- lapply(names(sets), function(pathway) {
    members <- sets[[pathway]]
    non_members <- setdiff(background, members)
    member_scores <- feature_scores[members]
    background_scores <- feature_scores[non_members]

    p <- NA_real_
    if (length(member_scores) > 0 && length(background_scores) > 0) {
      p <- stats::wilcox.test(
        member_scores,
        background_scores,
        alternative = alternative,
        exact = FALSE
      )$p.value
    }
    delta <- stats::median(member_scores, na.rm = TRUE) - stats::median(background_scores, na.rm = TRUE)

    data.frame(
      pathway = pathway,
      database = meta$database[match(pathway, meta$pathway)],
      n_background = length(background),
      n_pathway = length(members),
      n_non_pathway = length(non_members),
      statistic = mean(member_scores, na.rm = TRUE),
      median_pathway = stats::median(member_scores, na.rm = TRUE),
      median_background = stats::median(background_scores, na.rm = TRUE),
      effect = delta,
      dir = sign(delta),
      p = p,
      alternative = alternative,
      pathway_genes = paste(members, collapse = ";"),
      method = "wilcoxon",
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, res)
  if (isTRUE(fdr_by_database)) {
    out <- adjust_fdr_by_database(out, p_col = "p", database_col = "database", method = p_adjust_method)
    out$fdr_scope <- "database"
  } else {
    out$fdr <- stats::p.adjust(out$p, method = p_adjust_method)
    out$fdr_scope <- "global"
  }
  out$signed.p <- signed_log10_p(out$p, out$dir)
  out$signed.fdr <- signed_log10_p(out$fdr, out$dir)
  out[order(out$p, -abs(out$effect)), ]
}

#' Thresholded Wilcoxon Pathway Enrichment
#'
#' Score-thresholded Wilcoxon enrichment.
#' first restrict the scored feature universe to features with absolute score at
#' or above `threshold`, then run Wilcoxon enrichment on that thresholded
#' background. This keeps the statistical logic explicit rather than silently
#' changing the background.
#'
#' @param feature_scores Named numeric vector.
#' @param pathway_db An `EnrichKit_pathway_db` object or named list.
#' @param threshold Absolute score threshold.
#' @param background Optional analysis-specific background.
#' @param ... Additional arguments passed to [wilcox_enrichment()].
#'
#' @return A Wilcoxon enrichment table with threshold metadata.
#' @export
thresholded_wilcox_enrichment <- function(feature_scores,
                                          pathway_db,
                                          threshold,
                                          background = NULL,
                                          ...) {
  feature_scores <- clean_named_numeric(feature_scores)
  if (is.null(background)) {
    background <- names(feature_scores)
  }
  background <- intersect(unique(as.character(background)), names(feature_scores))
  threshold_background <- background[abs(feature_scores[background]) >= threshold]
  if (length(threshold_background) < 2) {
    stop("Thresholded background contains fewer than two scored features.")
  }
  out <- wilcox_enrichment(
    feature_scores = feature_scores,
    pathway_db = pathway_db,
    background = threshold_background,
    ...
  )
  out$threshold <- threshold
  out$n_threshold_background <- length(threshold_background)
  out$method <- "thresholded_wilcoxon"
  out
}

#' Rank-Shift Wilcoxon Pathway Enrichment
#'
#' Rank-shift Wilcoxon enrichment, a method from Xiaoyu Song. Scores are
#' converted to ranks, non-pathway ranks are shifted by a user-defined fraction
#' of the ranked feature universe, and one-sided Wilcoxon tests are performed
#' separately for enrichment and depletion. This tests whether pathway-member
#' ranks remain higher or lower than non-pathway ranks even after imposing a
#' rank-shift margin on the non-pathway background.
#'
#' @param feature_scores Named numeric vector of feature-level scores.
#' @param pathway_db An `EnrichKit_pathway_db` object or named list.
#' @param background Optional analysis-specific background features.
#' @param shift_fraction Fraction of scored background features used as the rank
#'   shift margin. Xiaoyu Song's example used `0.20`.
#' @param direction Direction to test. `"enriched"` runs a one-sided greater
#'   test against non-pathway ranks shifted upward; `"depleted"` runs a one-sided
#'   less test against non-pathway ranks shifted downward; `"both"` returns both.
#' @param min_size Minimum pathway size after background matching.
#' @param max_size Maximum pathway size after background matching.
#' @param p_adjust_method P-value adjustment method.
#' @param fdr_by_database If `TRUE`, adjust FDR separately within each database
#'   label. If `FALSE`, adjust across all returned pathways together.
#'
#' @return A data frame of pathway-level rank-shift Wilcoxon results.
#' @export
rank_shift_wilcox_enrichment <- function(feature_scores,
                                         pathway_db,
                                         background = NULL,
                                         shift_fraction = 0.20,
                                         direction = c("both", "enriched", "depleted"),
                                         min_size = 1,
                                         max_size = Inf,
                                         p_adjust_method = "BH",
                                         fdr_by_database = TRUE) {
  direction <- match.arg(direction)
  if (!is.numeric(shift_fraction) || length(shift_fraction) != 1 ||
      !is.finite(shift_fraction) || shift_fraction < 0) {
    stop("`shift_fraction` must be a single non-negative number.")
  }

  feature_scores <- clean_named_numeric(feature_scores)
  if (is.null(background)) {
    background <- names(feature_scores)
  }
  background <- intersect(unique(as.character(background)), names(feature_scores))
  if (length(background) < 2) {
    stop("Background must contain at least two scored features.")
  }

  score_rank <- rank(feature_scores[background], ties.method = "average", na.last = "keep")
  rank_shift <- round(length(background) * shift_fraction)
  sets <- clean_gene_sets(as_gene_sets(pathway_db), universe = background, min_size = min_size, max_size = max_size)
  if (length(sets) == 0) {
    return(data.frame())
  }
  meta <- get_pathway_meta_for_results(pathway_db, names(sets))

  directions <- if (direction == "both") c("enriched", "depleted") else direction
  res <- lapply(directions, function(test_direction) {
    do.call(rbind, lapply(names(sets), function(pathway) {
      members <- sets[[pathway]]
      non_members <- setdiff(background, members)
      score_in <- score_rank[members]
      score_out <- score_rank[non_members]

      if (test_direction == "enriched") {
        shifted_out <- score_out + rank_shift
        alternative <- "greater"
        dir <- 1
      } else {
        shifted_out <- score_out - rank_shift
        alternative <- "less"
        dir <- -1
      }

      p <- NA_real_
      if (length(score_in) > 0 && length(shifted_out) > 0) {
        p <- stats::wilcox.test(
          score_in,
          shifted_out,
          alternative = alternative,
          exact = FALSE
        )$p.value
      }

      data.frame(
        pathway = pathway,
        database = meta$database[match(pathway, meta$pathway)],
        n_background = length(background),
        n_pathway = length(members),
        n_non_pathway = length(non_members),
        shift_fraction = shift_fraction,
        rank_shift = rank_shift,
        direction_test = test_direction,
        alternative = alternative,
        median_pathway_rank = stats::median(score_in, na.rm = TRUE),
        median_shifted_background_rank = stats::median(shifted_out, na.rm = TRUE),
        effect = stats::median(score_in, na.rm = TRUE) - stats::median(shifted_out, na.rm = TRUE),
        dir = dir,
        p = p,
        pathway_genes = paste(members, collapse = ";"),
        method = "rank_shift_wilcoxon_xiaoyu_song",
        stringsAsFactors = FALSE
      )
    }))
  })

  out <- do.call(rbind, res)
  if (isTRUE(fdr_by_database)) {
    out$fdr <- NA_real_
    split_key <- interaction(out$database, out$direction_test, drop = TRUE)
    for (idx in split(seq_len(nrow(out)), split_key)) {
      out$fdr[idx] <- stats::p.adjust(out$p[idx], method = p_adjust_method)
    }
    out$fdr_scope <- "database_by_direction"
  } else {
    out$fdr <- stats::p.adjust(out$p, method = p_adjust_method)
    out$fdr_scope <- "global"
  }
  out$signed.p <- signed_log10_p(out$p, out$dir)
  out$signed.fdr <- signed_log10_p(out$fdr, out$dir)
  out[order(out$p, out$pathway), ]
}

#' Database-Specific FDR Adjustment
#'
#' Adjust p-values within each pathway database. This preserves the
#' database-wise correction pattern used in the HOPE pathway workflows.
#'
#' @param results Enrichment result data frame.
#' @param p_col P-value column.
#' @param database_col Database label column.
#' @param method Adjustment method passed to [stats::p.adjust()].
#'
#' @return `results` with `fdr` and `signed.fdr` columns updated/added.
#' @export
adjust_fdr_by_database <- function(results,
                                   p_col = "p",
                                   database_col = "database",
                                   method = "BH") {
  if (!p_col %in% colnames(results)) {
    stop("Missing p-value column: ", p_col)
  }
  if (!database_col %in% colnames(results)) {
    stop("FDR adjustment by database requested without database labels.")
  }
  split_idx <- split(seq_len(nrow(results)), results[[database_col]], drop = TRUE)
  fdr <- rep(NA_real_, nrow(results))
  for (idx in split_idx) {
    fdr[idx] <- stats::p.adjust(results[[p_col]][idx], method = method)
  }
  results$fdr <- fdr
  if ("dir" %in% colnames(results)) {
    results$signed.fdr <- signed_log10_p(results$fdr, results$dir)
  }
  results
}

#' Signed -log10 P-value
#'
#' @param p Numeric p-values.
#' @param direction Numeric direction/effect sign.
#'
#' @return Signed `-log10(p)`.
#' @export
signed_log10_p <- function(p, direction) {
  -log10(pmax(as.numeric(p), .Machine$double.xmin)) * sign(as.numeric(direction))
}

clean_named_numeric <- function(x) {
  if (is.data.frame(x)) {
    stop("Use `rank_features()` before passing data frames to enrichment functions.")
  }
  if (is.null(names(x))) {
    stop("Scores must be a named numeric vector.")
  }
  out <- as.numeric(x)
  names(out) <- names(x)
  out <- out[!is.na(names(out)) & nzchar(names(out)) & is.finite(out)]
  out
}

get_pathway_meta_for_results <- function(pathway_db, pathways) {
  if (inherits(pathway_db, "EnrichKit_pathway_db")) {
    meta <- pathway_db$metadata
    meta[match(pathways, meta$pathway), c("pathway", "database"), drop = FALSE]
  } else {
    data.frame(pathway = pathways, database = "custom", stringsAsFactors = FALSE)
  }
}
