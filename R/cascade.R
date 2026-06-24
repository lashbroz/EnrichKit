#' Build Cascade Membership Data
#'
#' Creates a long table showing which genes/features appear in each pathway.
#' This is useful for cascade-style plots that track shared members across a
#' ranked sequence of related pathways.
#'
#' @param gene_sets Named list of gene sets or an `EnrichKit_pathway_db` object.
#' @param pathways Optional pathway order/subset.
#' @param members Optional member order/subset.
#'
#' @return A data frame with `pathway`, `member`, `present`, and order columns.
#' @export
cascade_membership_data <- function(gene_sets,
                                     pathways = NULL,
                                     members = NULL) {
  sets <- as_gene_sets(gene_sets)
  if (is.null(pathways)) {
    pathways <- names(sets)
  }
  pathways <- intersect(as.character(pathways), names(sets))
  if (length(pathways) == 0) {
    return(data.frame())
  }
  if (is.null(members)) {
    members <- unique(unlist(sets[pathways], use.names = FALSE))
  }
  members <- as.character(members)

  grid <- expand.grid(
    pathway = pathways,
    member = members,
    stringsAsFactors = FALSE
  )
  grid$present <- mapply(function(pathway, member) member %in% sets[[pathway]], grid$pathway, grid$member)
  grid$pathway_order <- match(grid$pathway, pathways)
  grid$member_order <- match(grid$member, members)
  grid
}

#' Plot Cascade Membership
#'
#' Base-R cascade plot for pathway membership. The function returns the plotted
#' data invisibly so the same table can be passed to custom plotting code.
#'
#' @param cascade_data Output from [cascade_membership_data()] or a gene-set
#'   object accepted by [cascade_membership_data()].
#' @param main Plot title.
#' @param present_col Color for present memberships.
#' @param absent_col Color for absent memberships.
#' @param cex_axis Axis text size.
#'
#' @return Invisibly returns the cascade data frame.
#' @export
plot_cascade_membership <- function(cascade_data,
                                    main = "Pathway membership cascade",
                                    present_col = "#2b8cbe",
                                    absent_col = "grey92",
                                    cex_axis = 0.7) {
  if (!is.data.frame(cascade_data)) {
    cascade_data <- cascade_membership_data(cascade_data)
  }
  if (nrow(cascade_data) == 0) {
    stop("No cascade data to plot.")
  }
  pathways <- unique(cascade_data$pathway[order(cascade_data$pathway_order)])
  members <- unique(cascade_data$member[order(cascade_data$member_order)])
  mat <- matrix(FALSE, nrow = length(pathways), ncol = length(members), dimnames = list(pathways, members))
  idx <- cascade_data$present
  mat[cbind(cascade_data$pathway_order[idx], cascade_data$member_order[idx])] <- TRUE

  graphics::image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1, , drop = FALSE]),
    col = c(absent_col, present_col),
    axes = FALSE,
    xlab = "Member",
    ylab = "Pathway",
    main = main
  )
  graphics::axis(1, at = seq_along(members), labels = members, las = 2, cex.axis = cex_axis)
  graphics::axis(2, at = seq_along(pathways), labels = rev(pathways), las = 2, cex.axis = cex_axis)
  graphics::box()
  invisible(cascade_data)
}

#' Cascade Threshold Pathway Consolidation
#'
#' Consolidates an ordered pathway result table using a cascade threshold rule
#' developed in the CPTAC/CBTN workflow. Pathways are converted to a binary
#' pathway-by-member matrix, columns are cascade-ordered by shared membership,
#' and pathways are walked in priority order. The first pathway is retained;
#' later pathways are retained only if they introduce at least
#' `min_new_members` genes/features not already encountered in higher-priority
#' pathways. This is useful when many significant pathways contain nearly the
#' same genes but a few later pathways add a distinct block of biology.
#'
#' @param results Ordered enrichment result table.
#' @param gene_sets Named list of gene sets or an `EnrichKit_pathway_db` object.
#' @param pathway_col Column in `results` containing pathway names.
#' @param members Optional member universe to display/test. If `NULL`, all
#'   members in selected pathways are used.
#' @param min_new_members Minimum number of newly introduced members required
#'   to retain a non-first pathway.
#' @param row_order `"input"` preserves result order; `"size_desc"` orders
#'   pathways by matched set size before applying the cascade rule.
#'
#' @return A list with `results`, `membership_matrix`, and `retained_pathways`.
#' @export
cascade_threshold_consolidation <- function(results,
                                            gene_sets,
                                            pathway_col = "pathway",
                                            members = NULL,
                                            min_new_members = 5,
                                            row_order = c("input", "size_desc")) {
  row_order <- match.arg(row_order)
  if (!pathway_col %in% colnames(results)) {
    stop("Missing pathway column: ", pathway_col)
  }
  if (!is.numeric(min_new_members) || length(min_new_members) != 1 || min_new_members < 0) {
    stop("`min_new_members` must be a single non-negative number.")
  }

  sets <- as_gene_sets(gene_sets)
  pathways <- intersect(as.character(results[[pathway_col]]), names(sets))
  if (length(pathways) == 0) {
    stop("No result pathways were found in `gene_sets`.")
  }
  if (is.null(members)) {
    members <- unique(unlist(sets[pathways], use.names = FALSE))
  }
  members <- as.character(members)

  mat <- pathway_member_matrix(sets, pathways = pathways, members = members)
  mat <- mat[rowSums(mat) > 0, colSums(mat) > 0, drop = FALSE]
  if (row_order == "size_desc") {
    mat <- mat[order(rowSums(mat), decreasing = TRUE), , drop = FALSE]
  }
  mat <- cascade_order_membership_matrix(mat)

  seen <- character(0)
  new_members <- vector("list", nrow(mat))
  n_new <- integer(nrow(mat))
  retained <- logical(nrow(mat))
  for (i in seq_len(nrow(mat))) {
    present <- colnames(mat)[mat[i, ] > 0]
    new <- setdiff(present, seen)
    new_members[[i]] <- new
    n_new[i] <- length(new)
    retained[i] <- i == 1 || n_new[i] >= min_new_members
    seen <- union(seen, present)
  }

  out <- results[match(rownames(mat), results[[pathway_col]]), , drop = FALSE]
  out$cascade_order <- seq_len(nrow(mat))
  out$cascade_n_members <- rowSums(mat)
  out$cascade_n_new_members <- n_new
  out$cascade_new_members <- vapply(new_members, paste, character(1), collapse = ";")
  out$cascade_retained <- retained
  out$cascade_threshold <- min_new_members
  rownames(out) <- NULL

  list(
    results = out,
    membership_matrix = mat,
    retained_pathways = out[[pathway_col]][out$cascade_retained]
  )
}

#' Plot Cascade Threshold Consolidation
#'
#' Base-R visualization for [cascade_threshold_consolidation()]. Rows are
#' pathways, columns are genes/features, black cells indicate membership, and
#' retained pathways are marked on the left.
#'
#' @param cascade Output from [cascade_threshold_consolidation()].
#' @param main Plot title.
#' @param retained_col Color used for retained pathway markers.
#' @param member_col Color used for member cells.
#' @param absent_col Color used for absent cells.
#' @param cex_axis Axis text size.
#'
#' @return Invisibly returns `cascade`.
#' @export
plot_cascade_threshold_consolidation <- function(cascade,
                                                 main = "Cascade threshold pathway consolidation",
                                                 retained_col = "#d7301f",
                                                 member_col = "black",
                                                 absent_col = "grey95",
                                                 cex_axis = 0.7) {
  if (!is.list(cascade) || !all(c("results", "membership_matrix") %in% names(cascade))) {
    stop("`cascade` must be output from `cascade_threshold_consolidation()`.")
  }
  mat <- cascade$membership_matrix
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop("No membership matrix to plot.")
  }
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1, 12))
  on.exit(graphics::layout(1), add = TRUE)

  retained <- cascade$results$cascade_retained
  marker <- matrix(as.numeric(retained), ncol = 1)
  graphics::image(
    x = 1,
    y = seq_len(nrow(mat)),
    z = t(marker[nrow(marker):1, , drop = FALSE]),
    col = c(absent_col, retained_col),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "keep"
  )
  graphics::box()

  graphics::image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1, , drop = FALSE]),
    col = c(absent_col, member_col),
    axes = FALSE,
    xlab = "Member",
    ylab = "Pathway",
    main = main
  )
  graphics::axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = cex_axis)
  graphics::axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = cex_axis)
  graphics::box()
  invisible(cascade)
}

pathway_member_matrix <- function(sets, pathways, members) {
  mat <- matrix(0L, nrow = length(pathways), ncol = length(members), dimnames = list(pathways, members))
  for (pathway in pathways) {
    mat[pathway, members %in% sets[[pathway]]] <- 1L
  }
  mat
}

cascade_order_membership_matrix <- function(mat) {
  if (nrow(mat) < 1 || ncol(mat) < 1) {
    return(mat)
  }
  mat <- mat[order(rowSums(mat), decreasing = TRUE), , drop = FALSE]
  pattern <- rep("", ncol(mat))
  for (i in seq_len(nrow(mat))) {
    pattern <- paste0(pattern, mat[i, ])
    ord <- order(pattern, decreasing = TRUE)
    mat <- mat[, ord, drop = FALSE]
    pattern <- pattern[ord]
  }
  mat[, colSums(mat) > 0, drop = FALSE]
}
