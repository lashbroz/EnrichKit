#' Clean Gene Sets
#'
#' Clean a named list of feature or gene sets by removing missing/empty IDs,
#' deduplicating members within each set, optionally intersecting with a
#' background universe, and filtering by set size.
#'
#' @param gene_sets Named list of character vectors.
#' @param universe Optional character vector of valid analysis features.
#' @param min_size Minimum retained set size after universe filtering.
#' @param max_size Maximum retained set size after universe filtering.
#'
#' @return A named list of cleaned gene sets.
#' @export
clean_gene_sets <- function(gene_sets,
                            universe = NULL,
                            min_size = 1,
                            max_size = Inf) {
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("`gene_sets` must be a named list.")
  }
  if (any(is.na(names(gene_sets)) | !nzchar(names(gene_sets)))) {
    stop("All gene sets must have non-empty names.")
  }
  if (anyDuplicated(names(gene_sets))) {
    stop("Duplicate pathway names detected: ",
         paste(unique(names(gene_sets)[duplicated(names(gene_sets))]), collapse = ", "))
  }

  if (!is.null(universe)) {
    universe <- unique(as.character(universe))
    universe <- universe[!is.na(universe) & nzchar(universe)]
  }

  out <- lapply(gene_sets, function(x) {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    if (!is.null(universe)) {
      x <- intersect(x, universe)
    }
    sort(x)
  })

  set_size <- lengths(out)
  out[set_size >= min_size & set_size <= max_size]
}

#' Create a Pathway Database Object
#'
#' Create a structured EnrichKit pathway database from a user-provided named
#' list of gene sets. This is the main entry point for custom pathway lists,
#' MSigDB-derived lists, kinase substrate lists, glyco sets, SUMER modules, or
#' manuscript-specific curated signatures.
#'
#' @param gene_sets Named list mapping pathway names to gene/feature IDs.
#' @param database Either a scalar database label applied to all sets, or a
#'   named character vector with one label per pathway.
#' @param source Optional source label, such as `"MSigDB"` or `"manual"`.
#' @param version Optional source/database version.
#' @param universe Optional analysis background. If supplied, set members are
#'   intersected with the universe and dropped when outside size limits.
#' @param min_size Minimum retained set size.
#' @param max_size Maximum retained set size.
#' @param deduplicate_identical If `TRUE`, remove exactly duplicated gene sets,
#'   retaining the first occurrence.
#'
#' @return An object of class `EnrichKit_pathway_db`.
#' @export
make_pathway_db <- function(gene_sets,
                            database = "custom",
                            source = "user",
                            version = NA_character_,
                            universe = NULL,
                            min_size = 1,
                            max_size = Inf,
                            deduplicate_identical = TRUE) {
  original_n_sets <- length(gene_sets)
  cleaned <- clean_gene_sets(
    gene_sets = gene_sets,
    universe = universe,
    min_size = min_size,
    max_size = max_size
  )

  if (deduplicate_identical && length(cleaned) > 0) {
    collapsed <- vapply(cleaned, paste, character(1), collapse = "\r")
    cleaned <- cleaned[match(unique(collapsed), collapsed)]
  }

  if (length(cleaned) == 0) {
    warning("No gene sets remain after cleaning/filtering.")
  }

  database <- resolve_database_labels(database, names(cleaned))

  metadata <- data.frame(
    pathway = names(cleaned),
    database = unname(database[names(cleaned)]),
    source = source,
    version = version,
    n_genes = lengths(cleaned),
    stringsAsFactors = FALSE
  )

  out <- list(
    sets = cleaned,
    metadata = metadata,
    source = source,
    version = version,
    universe = universe,
    parameters = list(
      min_size = min_size,
      max_size = max_size,
      deduplicate_identical = deduplicate_identical,
      original_n_sets = original_n_sets,
      retained_n_sets = length(cleaned)
    )
  )
  class(out) <- c("EnrichKit_pathway_db", "list")
  validate_pathway_db(out)
  out
}

#' Extract Gene Sets from a Pathway Database
#'
#' @param pathway_db An `EnrichKit_pathway_db` object or a named list.
#'
#' @return A named list of gene sets.
#' @export
as_gene_sets <- function(pathway_db) {
  if (inherits(pathway_db, "EnrichKit_pathway_db")) {
    return(pathway_db$sets)
  }
  if (is.list(pathway_db) && !is.null(names(pathway_db))) {
    return(pathway_db)
  }
  stop("`pathway_db` must be an EnrichKit pathway database or named list.")
}

#' Extract Pathway Metadata
#'
#' @param pathway_db An `EnrichKit_pathway_db` object.
#'
#' @return A data frame with pathway metadata.
#' @export
pathway_metadata <- function(pathway_db) {
  validate_pathway_db(pathway_db)
  pathway_db$metadata
}

#' Match Gene Sets to an Analysis Background
#'
#' Intersect every pathway with an analysis-specific background and return a new
#' pathway database with updated gene counts. This matters because a small
#' interrogated/eligible feature universe can cause many pathways to collapse onto
#' the same few retained genes. Pathways are trimmed to the background, dropped
#' if their matched size is outside `min_size`/`max_size`, and the returned
#' object stores a `matching_summary` table so this transformation is auditable.
#' Members within each pathway are sorted alphabetically by [clean_gene_sets()].
#' Pathway order is source/input order by default, with optional deterministic
#' ordering via `order_by`.
#'
#' @param pathway_db An `EnrichKit_pathway_db` object or named list.
#' @param background Character vector of valid analysis features.
#' @param min_size Minimum retained set size after matching.
#' @param max_size Maximum retained set size after matching.
#' @param warn If `TRUE`, warn when any pathway loses members.
#' @param order_by Ordering for retained pathways. `"input"` preserves source
#'   order; `"pathway"` sorts alphabetically; `"database"` sorts by database
#'   then pathway; `"size_desc"`/`"size_asc"` sort by matched pathway size.
#'
#' @return An `EnrichKit_pathway_db` object matched to `background`.
#' @export
match_pathway_background <- function(pathway_db,
                                     background,
                                     min_size = 1,
                                     max_size = Inf,
                                     warn = TRUE,
                                     order_by = c("input", "pathway", "database", "size_desc", "size_asc")) {
  order_by <- match.arg(order_by)
  if (inherits(pathway_db, "EnrichKit_pathway_db")) {
    sets <- pathway_db$sets
    meta <- pathway_db$metadata
  } else {
    sets <- as_gene_sets(pathway_db)
    meta <- data.frame(
      pathway = names(sets),
      database = "custom",
      source = "user",
      version = NA_character_,
      n_genes = lengths(sets),
      stringsAsFactors = FALSE
    )
  }

  background <- unique(as.character(background))
  background <- background[!is.na(background) & nzchar(background)]
  original_counts <- lengths(sets)
  matched <- clean_gene_sets(sets, universe = background, min_size = min_size, max_size = max_size)
  matched_counts <- lengths(matched)
  matched_all <- lapply(sets, function(x) sort(unique(intersect(as.character(x), background))))
  matched_all_counts <- lengths(matched_all)
  summary <- data.frame(
    pathway = names(sets),
    database = meta$database[match(names(sets), meta$pathway)],
    original_n_genes = as.integer(original_counts),
    matched_n_genes = as.integer(matched_all_counts),
    n_genes_lost = as.integer(original_counts - matched_all_counts),
    retained = names(sets) %in% names(matched),
    drop_reason = ifelse(
      matched_all_counts < min_size,
      "below_min_size",
      ifelse(matched_all_counts > max_size, "above_max_size", "retained")
    ),
    stringsAsFactors = FALSE
  )

  if (warn) {
    common <- intersect(names(sets), names(matched))
    lost <- original_counts[common] - matched_counts[common]
    if (any(lost > 0)) {
      warning(sum(lost > 0), " pathway(s) lost genes after background matching.")
    }
    dropped <- setdiff(names(sets), names(matched))
    if (length(dropped) > 0) {
      warning(length(dropped), " pathway(s) dropped after background/size filtering.")
    }
  }

  meta <- meta[match(names(matched), meta$pathway), , drop = FALSE]
  meta$n_genes <- lengths(matched)
  ord <- order_matched_pathways(meta, order_by = order_by)
  meta <- meta[ord, , drop = FALSE]
  matched <- matched[meta$pathway]

  out <- list(
    sets = matched,
    metadata = meta,
    matching_summary = summary,
    source = if (nrow(meta) > 0) unique(meta$source)[1] else NA_character_,
    version = if (nrow(meta) > 0) unique(meta$version)[1] else NA_character_,
    universe = background,
    parameters = list(min_size = min_size, max_size = max_size, order_by = order_by)
  )
  class(out) <- c("EnrichKit_pathway_db", "list")
  validate_pathway_db(out)
  out
}

#' Extract Pathway Background-Matching Summary
#'
#' @param pathway_db An `EnrichKit_pathway_db` returned by
#'   [match_pathway_background()].
#'
#' @return A data frame with original counts, matched counts, retained status,
#'   and drop reasons.
#' @export
pathway_matching_summary <- function(pathway_db) {
  validate_pathway_db(pathway_db)
  if (is.null(pathway_db$matching_summary)) {
    stop("No matching summary found. Run `match_pathway_background()` first.")
  }
  pathway_db$matching_summary
}

#' Count Pathway Members in a Background
#'
#' @param pathway_db An `EnrichKit_pathway_db` object or named list.
#' @param background Optional character vector of valid analysis features.
#'
#' @return A data frame with pathway gene counts.
#' @export
pathway_gene_counts <- function(pathway_db, background = NULL) {
  sets <- as_gene_sets(pathway_db)
  if (!is.null(background)) {
    background <- unique(as.character(background))
  }
  data.frame(
    pathway = names(sets),
    n_genes = lengths(sets),
    n_in_background = if (is.null(background)) lengths(sets) else lengths(lapply(sets, intersect, background)),
    stringsAsFactors = FALSE
  )
}

resolve_database_labels <- function(database, pathways) {
  if (length(database) == 1 && is.null(names(database))) {
    out <- stats::setNames(rep(as.character(database), length(pathways)), pathways)
    return(out)
  }
  if (is.null(names(database))) {
    stop("When `database` is not scalar, it must be named by pathway.")
  }
  missing <- setdiff(pathways, names(database))
  if (length(missing) > 0) {
    stop("Missing database labels for pathways: ", paste(missing, collapse = ", "))
  }
  out <- as.character(database[pathways])
  names(out) <- pathways
  out
}

order_matched_pathways <- function(meta, order_by = "input") {
  if (nrow(meta) == 0) {
    return(integer(0))
  }
  if (order_by == "input") {
    return(seq_len(nrow(meta)))
  }
  if (order_by == "pathway") {
    return(order(meta$pathway))
  }
  if (order_by == "database") {
    return(order(meta$database, meta$pathway))
  }
  if (order_by == "size_desc") {
    return(order(-meta$n_genes, meta$pathway))
  }
  if (order_by == "size_asc") {
    return(order(meta$n_genes, meta$pathway))
  }
  seq_len(nrow(meta))
}

validate_pathway_db <- function(pathway_db) {
  if (!inherits(pathway_db, "EnrichKit_pathway_db")) {
    stop("Expected an EnrichKit pathway database object.")
  }
  if (!is.list(pathway_db$sets) || is.null(names(pathway_db$sets))) {
    stop("`pathway_db$sets` must be a named list.")
  }
  if (!is.data.frame(pathway_db$metadata)) {
    stop("`pathway_db$metadata` must be a data frame.")
  }
  required <- c("pathway", "database", "n_genes")
  missing <- setdiff(required, colnames(pathway_db$metadata))
  if (length(missing) > 0) {
    stop("Pathway metadata missing required columns: ", paste(missing, collapse = ", "))
  }
  if (!identical(names(pathway_db$sets), pathway_db$metadata$pathway)) {
    stop("Pathway set names must match metadata$pathway in order.")
  }
  invisible(TRUE)
}
