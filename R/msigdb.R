#' Read a GMT Gene-Set File
#'
#' Reads a standard GMT file into a named list of gene sets. Each row is
#' interpreted as pathway name, description, then members.
#'
#' @param file GMT file path.
#' @param min_size Minimum retained set size.
#' @param max_size Maximum retained set size.
#' @param universe Optional analysis universe used to filter members.
#'
#' @return A named list of character vectors.
#' @export
read_gmt <- function(file,
                     min_size = 1,
                     max_size = Inf,
                     universe = NULL) {
  if (!file.exists(file)) {
    stop("GMT file does not exist: ", file)
  }
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(lines)]
  sets <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3) {
      return(character(0))
    }
    fields[-c(1, 2)]
  })
  names(sets) <- vapply(strsplit(lines, "\t", fixed = TRUE), `[`, character(1), 1)
  clean_gene_sets(sets, universe = universe, min_size = min_size, max_size = max_size)
}

#' Build a Pathway Database from GMT
#'
#' @param file GMT file path.
#' @param database Database label to assign to all GMT pathways.
#' @param source Source label.
#' @param version Optional version label.
#' @param universe Optional analysis universe.
#' @param min_size Minimum retained set size.
#' @param max_size Maximum retained set size.
#'
#' @return An `EnrichKit_pathway_db` object.
#' @export
make_pathway_db_from_gmt <- function(file,
                                     database = tools::file_path_sans_ext(basename(file)),
                                     source = "GMT",
                                     version = NA_character_,
                                     universe = NULL,
                                     min_size = 1,
                                     max_size = Inf) {
  sets <- read_gmt(file, min_size = min_size, max_size = max_size, universe = universe)
  make_pathway_db(
    sets,
    database = database,
    source = source,
    version = version,
    universe = universe,
    min_size = min_size,
    max_size = max_size
  )
}

#' MSigDB Collection Metadata
#'
#' Returns common MSigDB GMT file names and collection labels. The function is
#' intentionally metadata-only; downloading requires explicit URLs or local
#' files because MSigDB access rules can vary by release.
#'
#' @param species One of `"human"` or `"mouse"`.
#' @param version MSigDB version label used to construct expected GMT names.
#'
#' @return A data frame with collection labels and GMT file names.
#' @export
msigdb_collections <- function(species = c("human", "mouse"),
                               version = NULL) {
  species <- match.arg(species)
  if (species == "human") {
    if (is.null(version)) {
      version <- "v2026.1.Hs"
    }
    out <- data.frame(
      collection = c("H", paste0("C", 1:9)),
      label = c(
        "Hallmark",
        "Positional",
        "Curated",
        "Regulatory target",
        "Computational",
        "Ontology",
        "Oncogenic signature",
        "Immunologic signature",
        "Cell type signature",
        "Immune signature"
      ),
      file_prefix = c("h", paste0("c", 1:9)),
      stringsAsFactors = FALSE
    )
  } else {
    if (is.null(version)) {
      version <- "v2026.1.Mm"
    }
    out <- data.frame(
      collection = c("MH", "M1", "M2", "M3", "M5", "M7", "M8"),
      label = c(
        "Mouse Hallmark",
        "Mouse positional",
        "Mouse curated",
        "Mouse regulatory target",
        "Mouse ontology",
        "Mouse immunologic signature",
        "Mouse cell type signature"
      ),
      file_prefix = c("mh", "m1", "m2", "m3", "m5", "m7", "m8"),
      stringsAsFactors = FALSE
    )
  }
  out$species <- species
  out$version <- version
  out$gmt_file <- paste0(out$file_prefix, ".all.", version, ".symbols.gmt")
  out[, c("collection", "label", "species", "version", "gmt_file")]
}

#' Expected MSigDB GMT File Names
#'
#' Returns the expected local GMT file names for selected MSigDB collections.
#' Users can use this as a checklist after downloading files from MSigDB.
#'
#' @param collections MSigDB collection codes. If `NULL`, returns all supported
#'   collections for `species`.
#' @param species One of `"human"` or `"mouse"`.
#' @param version Optional MSigDB version label.
#'
#' @return A named character vector of expected GMT file names.
#' @export
msigdb_expected_files <- function(collections = NULL,
                                  species = c("human", "mouse"),
                                  version = NULL) {
  tab <- msigdb_collections(species = species, version = version)
  if (is.null(collections)) {
    collections <- tab$collection
  }
  collections <- as.character(collections)
  missing <- setdiff(collections, tab$collection)
  if (length(missing) > 0) {
    stop("Unknown MSigDB collection(s): ", paste(missing, collapse = ", "))
  }
  out <- tab$gmt_file[match(collections, tab$collection)]
  names(out) <- collections
  out
}

#' Build a Combined MSigDB Pathway Database
#'
#' Reads one or more local MSigDB GMT files and combines them into a single
#' `EnrichKit_pathway_db`. Each collection is kept as its own database label so
#' downstream FDR correction can be database-wise.
#'
#' @param files Named character vector or named list of GMT file paths. Names
#'   should be MSigDB collection codes such as `"H"` or `"C2"`. If unnamed,
#'   collection labels are inferred from file names when possible.
#' @param species One of `"human"` or `"mouse"`.
#' @param version Optional MSigDB version label.
#' @param universe Optional analysis universe.
#' @param min_size Minimum retained set size.
#' @param max_size Maximum retained set size.
#' @param deduplicate_identical Remove exactly duplicated gene sets.
#'
#' @return An `EnrichKit_pathway_db` object.
#' @export
build_msigdb_pathway_db <- function(files,
                                    species = c("human", "mouse"),
                                    version = NULL,
                                    universe = NULL,
                                    min_size = 1,
                                    max_size = Inf,
                                    deduplicate_identical = TRUE) {
  species <- match.arg(species)
  files <- normalize_named_gmt_files(files, species = species, version = version)
  collections <- names(files)
  tab <- msigdb_collections(species = species, version = version)

  sets_by_collection <- lapply(files, read_gmt, min_size = min_size, max_size = max_size, universe = universe)
  all_sets <- unlist(sets_by_collection, recursive = FALSE, use.names = FALSE)
  collection_names <- rep(collections, lengths(sets_by_collection))
  names(all_sets) <- make.unique(unlist(lapply(sets_by_collection, names), use.names = FALSE), sep = "__dup")

  database <- paste0(
    "MSigDB_",
    collection_names,
    "_",
    tab$version[match(collection_names, tab$collection)]
  )
  names(database) <- names(all_sets)

  make_pathway_db(
    all_sets,
    database = database,
    source = "MSigDB",
    version = unique(tab$version[match(collection_names, tab$collection)]),
    universe = universe,
    min_size = min_size,
    max_size = max_size,
    deduplicate_identical = deduplicate_identical
  )
}

#' Build the KidsFirst Default Pathway Database
#'
#' Recreates the current KidsFirst pathway-database convention: use the HOPE
#' pathway GMT with `KEGG_MEDICUS*` pathways removed, add canonical MSigDB KEGG
#' pathways not already present, and filter to the provided KidsFirst/measured
#' gene universe.
#'
#' @param hope_gmt HOPE pathway database GMT.
#' @param canonical_kegg_gmt Canonical MSigDB KEGG GMT.
#' @param universe Measured gene universe used for filtering.
#' @param min_size Minimum retained set size. Defaults to KidsFirst usage, 6.
#' @param max_size Maximum retained set size. Defaults to KidsFirst usage, 249.
#'
#' @return An `EnrichKit_pathway_db` object.
#' @export
build_kfirst_default_pathway_db <- function(hope_gmt,
                                            canonical_kegg_gmt,
                                            universe,
                                            min_size = 6,
                                            max_size = 249) {
  hope_sets <- read_gmt(hope_gmt, min_size = 1, max_size = Inf)
  kegg_sets <- read_gmt(canonical_kegg_gmt, min_size = 1, max_size = Inf)

  hope_no_medicus <- hope_sets[!grepl("^KEGG_MEDICUS", names(hope_sets))]
  canonical_to_add <- kegg_sets[!names(kegg_sets) %in% names(hope_no_medicus)]
  combined <- c(hope_no_medicus, canonical_to_add)

  database <- c(
    stats::setNames(rep("HOPE_pathway_database_without_KEGG_MEDICUS", length(hope_no_medicus)), names(hope_no_medicus)),
    stats::setNames(rep("MSigDB_c2_cp_kegg_canonical", length(canonical_to_add)), names(canonical_to_add))
  )

  make_pathway_db(
    combined,
    database = database,
    source = "KidsFirst_default",
    version = NA_character_,
    universe = universe,
    min_size = min_size,
    max_size = max_size,
    deduplicate_identical = FALSE
  )
}

#' Load an Existing KidsFirst Gene-Set GMT
#'
#' Convenience wrapper for the already-built `gosets_all_kfirst.gmt` style file.
#'
#' @param file Existing KidsFirst GMT file.
#' @param source_table Optional source table with pathway/source columns.
#' @param universe Optional interrogated gene universe.
#' @param min_size Minimum retained set size.
#' @param max_size Maximum retained set size.
#'
#' @return An `EnrichKit_pathway_db` object.
#' @export
load_kfirst_gosets_gmt <- function(file,
                                   source_table = NULL,
                                   universe = NULL,
                                   min_size = 6,
                                   max_size = 249) {
  sets <- read_gmt(file, min_size = min_size, max_size = max_size, universe = universe)
  database <- stats::setNames(rep("KidsFirst_gosets_all", length(sets)), names(sets))
  if (!is.null(source_table)) {
    source_df <- utils::read.delim(source_table, stringsAsFactors = FALSE)
    source_col <- intersect(c("source", "database"), colnames(source_df))[1]
    if (!all(c("pathway", source_col) %in% colnames(source_df))) {
      stop("`source_table` must contain `pathway` and `source`/`database` columns.")
    }
    idx <- match(names(sets), source_df$pathway)
    database <- stats::setNames(as.character(source_df[[source_col]][idx]), names(sets))
    database[is.na(database)] <- "KidsFirst_gosets_all"
  }
  make_pathway_db(
    sets,
    database = database,
    source = "KidsFirst",
    version = NA_character_,
    universe = universe,
    min_size = min_size,
    max_size = max_size,
    deduplicate_identical = FALSE
  )
}

#' KidsFirst Default Pathway Database Components
#'
#' Describes the pathway sources used by the current KidsFirst default gene-set
#' database. If a source table is supplied, counts are read from that file.
#'
#' @param source_table Optional KidsFirst source TSV with `pathway`, `source`,
#'   and `n_genes` columns.
#'
#' @return A data frame describing pathway sources.
#' @export
kfirst_pathway_database_components <- function(source_table = NULL) {
  if (!is.null(source_table)) {
    x <- utils::read.delim(source_table, stringsAsFactors = FALSE)
    source_col <- intersect(c("source", "database"), colnames(x))[1]
    if (!all(c("pathway", source_col) %in% colnames(x))) {
      stop("`source_table` must contain `pathway` and `source`/`database` columns.")
    }
    counts <- as.data.frame(table(x[[source_col]]), stringsAsFactors = FALSE)
    colnames(counts) <- c("source", "n_pathways")
    counts$description <- kfirst_source_description(counts$source)
    return(counts[, c("source", "description", "n_pathways")])
  }
  data.frame(
    source = c("HOPE_pathway_database_without_KEGG_MEDICUS", "MSigDB_c2_cp_kegg_v7_canonical"),
    description = c(
      "HOPE pathway database after removing KEGG_MEDICUS* pathways",
      "Canonical KEGG pathways from MSigDB C2 canonical pathways"
    ),
    n_pathways = c(NA_integer_, NA_integer_),
    stringsAsFactors = FALSE
  )
}

kfirst_source_description <- function(source) {
  out <- rep("KidsFirst pathway source", length(source))
  out[source == "HOPE_pathway_database_without_KEGG_MEDICUS"] <-
    "HOPE pathway database after removing KEGG_MEDICUS* pathways"
  out[source %in% c("MSigDB_c2_cp_kegg_v7_canonical", "MSigDB_c2_cp_kegg_canonical")] <-
    "Canonical KEGG pathways from MSigDB C2 canonical pathways"
  out
}

#' Download a Gene-Set File
#'
#' Small explicit wrapper around [utils::download.file()]. This is useful for
#' Broad/MSigDB or other public GMT files when the URL is known.
#'
#' @param url Source URL.
#' @param destfile Destination file path.
#' @param overwrite Overwrite an existing destination.
#' @param mode Download mode passed to [utils::download.file()].
#'
#' @return Invisibly returns `destfile`.
#' @export
download_gene_set_file <- function(url,
                                   destfile,
                                   overwrite = FALSE,
                                   mode = "wb") {
  if (file.exists(destfile) && !overwrite) {
    stop("Destination already exists: ", destfile)
  }
  utils::download.file(url, destfile = destfile, mode = mode)
  invisible(destfile)
}

normalize_named_gmt_files <- function(files,
                                      species = c("human", "mouse"),
                                      version = NULL) {
  species <- match.arg(species)
  if (is.list(files)) {
    files <- unlist(files, use.names = TRUE)
  }
  files <- as.character(files)
  if (is.null(names(files)) || any(!nzchar(names(files)))) {
    inferred <- infer_msigdb_collection_from_file(files, species = species, version = version)
    if (any(is.na(inferred))) {
      stop("`files` must be named by MSigDB collection when collection cannot be inferred from file names.")
    }
    names(files) <- inferred
  }
  if (any(!file.exists(files))) {
    stop("Missing GMT file(s): ", paste(files[!file.exists(files)], collapse = ", "))
  }
  tab <- msigdb_collections(species = species, version = version)
  missing <- setdiff(names(files), tab$collection)
  if (length(missing) > 0) {
    stop("Unknown MSigDB collection name(s): ", paste(missing, collapse = ", "))
  }
  files
}

infer_msigdb_collection_from_file <- function(files,
                                              species = c("human", "mouse"),
                                              version = NULL) {
  tab <- msigdb_collections(species = species, version = version)
  lower_files <- tolower(basename(files))
  inferred <- rep(NA_character_, length(files))
  for (i in seq_along(files)) {
    hit <- tab$collection[tolower(tab$gmt_file) == lower_files[i]]
    if (length(hit) == 1) {
      inferred[i] <- hit
    }
  }
  inferred
}
