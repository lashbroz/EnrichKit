#' Prepare Enrichment Results for SUMER
#'
#' Writes the two standard SUMER inputs used in the HOPE/CPTAC workflow:
#' a pathway GMT file and a two-column pathway/weight file. This preserves the
#' usual behavior from `get_sumer.data`: if `weight_col` is absent, weights are
#' calculated as `-log10(fdr) * direction`.
#'
#' @param enrichment A data frame with pathway-level enrichment results.
#' @param gene_sets Named list of pathway gene sets.
#' @param out_prefix Output prefix, without extension.
#' @param pathway_col Column containing pathway names.
#' @param weight_col Optional column containing SUMER weights.
#' @param fdr_col FDR column used to derive weights when `weight_col` is absent.
#' @param direction_col Direction column used to derive weights when `weight_col`
#'   is absent. Values should be signed numeric direction/effect indicators.
#' @param min_abs_weight Optional minimum absolute weight to retain.
#' @param deduplicate Keep only the strongest row per pathway.
#'
#' @return A list containing written file paths and the filtered SUMER input.
#' @export
prepare_sumer_input <- function(enrichment,
                                gene_sets,
                                out_prefix,
                                pathway_col = "pathway",
                                weight_col = "weights",
                                fdr_col = "fdr",
                                direction_col = "dir",
                                min_abs_weight = NULL,
                                deduplicate = TRUE) {
  if (!is.data.frame(enrichment)) {
    stop("`enrichment` must be a data frame.")
  }
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("`gene_sets` must be a named list.")
  }
  if (!pathway_col %in% colnames(enrichment)) {
    stop("Missing pathway column: ", pathway_col)
  }

  x <- enrichment

  if (!weight_col %in% colnames(x)) {
    if (!all(c(fdr_col, direction_col) %in% colnames(x))) {
      stop("To derive weights, enrichment must contain `fdr_col` and `direction_col`.")
    }
    x[[weight_col]] <- -log10(pmax(as.numeric(x[[fdr_col]]), .Machine$double.xmin)) *
      sign(as.numeric(x[[direction_col]]))
  }

  sumer_data <- data.frame(
    pathway = as.character(x[[pathway_col]]),
    weights = as.numeric(x[[weight_col]]),
    stringsAsFactors = FALSE
  )
  sumer_data <- sumer_data[is.finite(sumer_data$weights), , drop = FALSE]
  sumer_data <- sumer_data[sumer_data$pathway %in% names(gene_sets), , drop = FALSE]

  if (!is.null(min_abs_weight)) {
    sumer_data <- sumer_data[abs(sumer_data$weights) >= min_abs_weight, , drop = FALSE]
  }

  sumer_data <- sumer_data[order(abs(sumer_data$weights), decreasing = TRUE), , drop = FALSE]
  if (deduplicate) {
    sumer_data <- sumer_data[match(unique(sumer_data$pathway), sumer_data$pathway), , drop = FALSE]
  }

  gmt_sets <- gene_sets[sumer_data$pathway]
  gmt_file <- paste0(out_prefix, "_pathways.gmt")
  data_file <- paste0(out_prefix, "_data.txt")

  write_gmt(gmt_sets, gmt_file)
  utils::write.table(
    sumer_data,
    file = data_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )

  list(
    gmt_file = gmt_file,
    data_file = data_file,
    sumer_data = sumer_data,
    gene_sets = gmt_sets
  )
}

#' Compatibility Wrapper for get_sumer.data
#'
#' Compatibility wrapper around [prepare_sumer_input()]. This preserves the
#' older HOPE/KidsFirst-style function name `get_sumer.data()` while routing the
#' implementation through the explicit EnrichKit SUMER input writer.
#'
#' Use this when migrating older analysis scripts that called `get_sumer.data`.
#' New code can call either `prepare_sumer_input()` or the snake-case alias
#' `get_sumer_data()`.
#'
#' @inheritParams prepare_sumer_input
#'
#' @return A list containing written file paths and the filtered SUMER input.
#' @export
get_sumer.data <- function(enrichment,
                           gene_sets,
                           out_prefix,
                           pathway_col = "pathway",
                           weight_col = "weights",
                           fdr_col = "fdr",
                           direction_col = "dir",
                           min_abs_weight = NULL,
                           deduplicate = TRUE) {
  prepare_sumer_input(
    enrichment = enrichment,
    gene_sets = gene_sets,
    out_prefix = out_prefix,
    pathway_col = pathway_col,
    weight_col = weight_col,
    fdr_col = fdr_col,
    direction_col = direction_col,
    min_abs_weight = min_abs_weight,
    deduplicate = deduplicate
  )
}

#' @rdname get_sumer.data
#' @export
get_sumer_data <- get_sumer.data

#' Write a SUMER Configuration Template
#'
#' Writes the SUMER configuration style used in the HOPE/CPTAC workflows:
#' `project`, `top_num`, and a `data` array containing `platform_name`,
#' `platform_abbr`, `gmt_file`, and `score_file`. SUMER is external and local
#' installations may support additional options, so the written file is intended
#' to be explicit and editable.
#'
#' @param sumer_input Output from [prepare_sumer_input()].
#' @param file Output config file path.
#' @param project SUMER project/output label.
#' @param top_num Number of top pathways SUMER should use per platform.
#' @param platform_name Human-readable platform name.
#' @param platform_abbr Short platform abbreviation used in SUMER node labels.
#' @param extra Named list of additional config fields to include.
#'
#' @return Invisibly returns `file`.
#' @export
write_sumer_config_template <- function(sumer_input,
                                        file,
                                        project = tools::file_path_sans_ext(basename(file)),
                                        top_num = 100,
                                        platform_name = project,
                                        platform_abbr = project,
                                        extra = list()) {
  required <- c("gmt_file", "data_file")
  missing <- setdiff(required, names(sumer_input))
  if (length(missing) > 0) {
    stop("`sumer_input` is missing required field(s): ", paste(missing, collapse = ", "))
  }
  if (!is.list(extra) || is.null(names(extra)) && length(extra) > 0) {
    stop("`extra` must be a named list.")
  }

  config <- c(
    list(
      project = project,
      top_num = top_num,
      data = list(list(
        platform_name = platform_name,
        platform_abbr = platform_abbr,
        gmt_file = normalizePath(sumer_input$gmt_file, mustWork = FALSE),
        score_file = normalizePath(sumer_input$data_file, mustWork = FALSE)
      ))
    ),
    extra
  )
  write_json_like_config(config, file)
  invisible(file)
}

#' Write a Named Gene-Set List as GMT
#'
#' @param gene_sets Named list of character vectors.
#' @param file Output GMT path.
#'
#' @return Invisibly returns `file`.
#' @export
write_gmt <- function(gene_sets, file) {
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("`gene_sets` must be a named list.")
  }
  lines <- vapply(names(gene_sets), function(nm) {
    genes <- unique(as.character(gene_sets[[nm]]))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    paste(c(nm, "na", genes), collapse = "\t")
  }, character(1))
  writeLines(lines, con = file)
  invisible(file)
}

#' Read SUMER Modules
#'
#' Reads SUMER `ap_sumer_edgelist.txt` and `ap_sumer_nodelist.txt` outputs and
#' returns connected-component modules with pathway weights, FDR, direction, set
#' size, degree, and betweenness when available. If no edge list exists, each
#' node is treated as its own module, matching the defensive behavior used in
#' the HOPE scripts.
#'
#' @param edge_file SUMER edge-list file, typically `ap_sumer_edgelist.txt`.
#' @param node_file SUMER node-list file, typically `ap_sumer_nodelist.txt`.
#' @param gene_sets Named list of pathway gene sets.
#' @param enrichment Optional enrichment table used to annotate nodes.
#' @param pathway_col Enrichment pathway column.
#' @param fdr_col Enrichment FDR column.
#' @param direction_col Enrichment direction column.
#' @param weight_col Enrichment/SUMER weight column.
#' @param strip_prefix Regex used to remove SUMER node prefixes.
#'
#' @return A list with `module_table`, `module_summary`, `module_list`, and
#'   `graph`.
#' @export
read_sumer_modules <- function(edge_file,
                               node_file,
                               gene_sets,
                               enrichment = NULL,
                               pathway_col = "pathway",
                               fdr_col = "fdr",
                               direction_col = "dir",
                               weight_col = "weights",
                               strip_prefix = "^.*?_([A-Z0-9])") {
  if (!file.exists(node_file)) {
    stop("Missing SUMER node file: ", node_file)
  }
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("`gene_sets` must be a named list.")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package `igraph` is required to read SUMER modules.")
  }

  node_data <- utils::read.table(node_file, header = TRUE, stringsAsFactors = FALSE)
  if (!"name" %in% colnames(node_data)) {
    stop("SUMER node file must contain a `name` column.")
  }
  node_names <- clean_sumer_node_names(node_data$name, strip_prefix = strip_prefix)

  edge_data <- tryCatch(
    utils::read.table(edge_file, header = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL,
    warning = function(w) NULL
  )

  if (!is.null(edge_data) && nrow(edge_data) > 0) {
    graph <- igraph::graph_from_data_frame(edge_data, directed = FALSE, vertices = node_data)
    igraph::V(graph)$name <- clean_sumer_node_names(igraph::V(graph)$name, strip_prefix = strip_prefix)
    module_list <- igraph::groups(igraph::components(graph))
    degree <- igraph::degree(graph, v = igraph::V(graph))
    between <- igraph::betweenness(graph, v = igraph::V(graph))
  } else {
    graph <- NULL
    module_list <- stats::setNames(as.list(node_names), seq_along(node_names))
    degree <- stats::setNames(rep(NA_real_, length(node_names)), node_names)
    between <- stats::setNames(rep(NA_real_, length(node_names)), node_names)
  }

  pathway_length <- lengths(gene_sets)
  annot <- make_sumer_node_annotation(
    pathways = unique(unlist(module_list, use.names = FALSE)),
    enrichment = enrichment,
    pathway_col = pathway_col,
    fdr_col = fdr_col,
    direction_col = direction_col,
    weight_col = weight_col
  )

  module_table <- do.call(rbind, lapply(seq_along(module_list), function(i) {
    nodes <- module_list[[i]]
    data.frame(
      module = i,
      pathway = nodes,
      n_genes = unname(pathway_length[match(nodes, names(pathway_length))]),
      degree = unname(degree[match(nodes, names(degree))]),
      betweenness = unname(between[match(nodes, names(between))]),
      stringsAsFactors = FALSE
    )
  }))

  module_table <- merge(module_table, annot, by = "pathway", all.x = TRUE, sort = FALSE)
  module_table <- module_table[order(module_table$module, module_table$fdr, -abs(module_table$weights)), ]
  rownames(module_table) <- NULL
  module_summary <- summarize_sumer_module_table(module_table)

  list(
    module_table = module_table,
    module_summary = module_summary,
    module_list = module_list,
    graph = graph
  )
}

#' Run SUMER
#'
#' Thin wrapper around the external SUMER function. SUMER is not reimplemented
#' by EnrichKit; this function only calls an installed/loaded SUMER function and
#' returns the expected output paths.
#'
#' @param config_file SUMER JSON configuration file.
#' @param output_name Output name/directory passed to SUMER.
#' @param sumer_fun Optional function to call. If omitted, EnrichKit searches
#'   for a loaded `sumer()` function.
#' @param expected_output_dir Directory expected to contain SUMER output files.
#' @param overwrite If `FALSE`, error when `expected_output_dir` already exists.
#'
#' @return Invisibly returns expected SUMER output paths.
#' @export
run_sumer <- function(config_file,
                      output_name,
                      sumer_fun = NULL,
                      expected_output_dir = output_name,
                      overwrite = FALSE) {
  if (!file.exists(config_file)) {
    stop("Missing SUMER config file: ", config_file)
  }

  if (dir.exists(expected_output_dir) && !overwrite) {
    stop(
      "SUMER output directory already exists: ", expected_output_dir,
      "\nUse overwrite = TRUE if you intend SUMER to overwrite/reuse it."
    )
  }

  if (is.null(sumer_fun)) {
    if (exists("sumer", mode = "function", inherits = TRUE)) {
      sumer_fun <- get("sumer", mode = "function", inherits = TRUE)
    } else {
      stop(
        "Could not find a loaded SUMER function. Load SUMER, or pass ",
        "`sumer_fun = sumer` explicitly.\n",
        "SUMER: https://github.com/bzhanglab/sumer"
      )
    }
  }

  sumer_fun(config_file, output_name)

  invisible(list(
    config_file = config_file,
    output_name = output_name,
    expected_output_dir = expected_output_dir,
    edge_file = file.path(expected_output_dir, "ap_sumer_edgelist.txt"),
    node_file = file.path(expected_output_dir, "ap_sumer_nodelist.txt")
  ))
}

#' Prepare, Run, and Optionally Read SUMER Results
#'
#' Convenience wrapper for the common workflow: write SUMER inputs from
#' enrichment results, run SUMER, and read the resulting modules.
#'
#' @inheritParams prepare_sumer_input
#' @inheritParams run_sumer
#' @param read_modules If `TRUE`, read SUMER modules after the run.
#'
#' @return A list with `prep`, `run`, and `modules`.
#' @export
prepare_run_read_sumer <- function(enrichment,
                                   gene_sets,
                                   out_prefix,
                                   config_file,
                                   output_name,
                                   pathway_col = "pathway",
                                   weight_col = "weights",
                                   fdr_col = "fdr",
                                   direction_col = "dir",
                                   min_abs_weight = NULL,
                                   deduplicate = TRUE,
                                   sumer_fun = NULL,
                                   expected_output_dir = output_name,
                                   overwrite = FALSE,
                                   read_modules = TRUE) {
  prep <- prepare_sumer_input(
    enrichment = enrichment,
    gene_sets = gene_sets,
    out_prefix = out_prefix,
    pathway_col = pathway_col,
    weight_col = weight_col,
    fdr_col = fdr_col,
    direction_col = direction_col,
    min_abs_weight = min_abs_weight,
    deduplicate = deduplicate
  )

  run <- run_sumer(
    config_file = config_file,
    output_name = output_name,
    sumer_fun = sumer_fun,
    expected_output_dir = expected_output_dir,
    overwrite = overwrite
  )

  modules <- NULL
  if (isTRUE(read_modules)) {
    if (!file.exists(run$node_file)) {
      warning("SUMER node file not found after run: ", run$node_file)
    } else {
      modules <- read_sumer_modules(
        edge_file = run$edge_file,
        node_file = run$node_file,
        gene_sets = gene_sets,
        enrichment = enrichment,
        pathway_col = pathway_col,
        fdr_col = fdr_col,
        direction_col = direction_col,
        weight_col = weight_col
      )
    }
  }

  list(prep = prep, run = run, modules = modules)
}

#' Prepare SUMER Inputs, Write Config Template, Run, and Read Modules
#'
#' High-level convenience wrapper for the full SUMER handoff. It writes SUMER
#' input files, writes a config template if requested, calls external SUMER, and
#' reads module outputs. Use `run = FALSE` to stop after creating the files.
#'
#' @inheritParams prepare_sumer_input
#' @param config_file Output or existing SUMER config file.
#' @param output_name Output name/directory passed to SUMER.
#' @param write_config If `TRUE`, write a config template before running.
#' @param config_extra Named list of additional config fields.
#' @param top_num Number of top pathways in the SUMER config template.
#' @param platform_name Human-readable platform name in the SUMER config.
#' @param platform_abbr Short platform abbreviation in the SUMER config.
#' @param run If `TRUE`, call [run_sumer()]. If `FALSE`, only prepare files.
#' @inheritParams run_sumer
#' @param read_modules If `TRUE`, read SUMER modules after running.
#'
#' @return A list with `prep`, `config_file`, `run`, and `modules`.
#' @export
sumer_workflow <- function(enrichment,
                           gene_sets,
                           out_prefix,
                           config_file = paste0(out_prefix, "_sumer_config.json"),
                           output_name = paste0(out_prefix, "_sumer_output"),
                           pathway_col = "pathway",
                           weight_col = "weights",
                           fdr_col = "fdr",
                           direction_col = "dir",
                           min_abs_weight = NULL,
                           deduplicate = TRUE,
                           write_config = TRUE,
                           config_extra = list(),
                           top_num = 100,
                           platform_name = output_name,
                           platform_abbr = output_name,
                           run = FALSE,
                           sumer_fun = NULL,
                           expected_output_dir = output_name,
                           overwrite = FALSE,
                           read_modules = TRUE) {
  prep <- prepare_sumer_input(
    enrichment = enrichment,
    gene_sets = gene_sets,
    out_prefix = out_prefix,
    pathway_col = pathway_col,
    weight_col = weight_col,
    fdr_col = fdr_col,
    direction_col = direction_col,
    min_abs_weight = min_abs_weight,
    deduplicate = deduplicate
  )

  if (isTRUE(write_config)) {
    write_sumer_config_template(
      sumer_input = prep,
      file = config_file,
      project = output_name,
      top_num = top_num,
      platform_name = platform_name,
      platform_abbr = platform_abbr,
      extra = config_extra
    )
  }

  run_out <- NULL
  modules <- NULL
  if (isTRUE(run)) {
    run_out <- run_sumer(
      config_file = config_file,
      output_name = output_name,
      sumer_fun = sumer_fun,
      expected_output_dir = expected_output_dir,
      overwrite = overwrite
    )
    if (isTRUE(read_modules)) {
      if (!file.exists(run_out$node_file)) {
        warning("SUMER node file not found after run: ", run_out$node_file)
      } else {
        modules <- read_sumer_modules(
          edge_file = run_out$edge_file,
          node_file = run_out$node_file,
          gene_sets = gene_sets,
          enrichment = enrichment,
          pathway_col = pathway_col,
          fdr_col = fdr_col,
          direction_col = direction_col,
          weight_col = weight_col
        )
      }
    }
  }

  list(
    prep = prep,
    config_file = config_file,
    run = run_out,
    modules = modules
  )
}

clean_sumer_node_names <- function(x, strip_prefix = "^.*?_([A-Z0-9])") {
  sub(strip_prefix, "\\1", as.character(x))
}

make_sumer_node_annotation <- function(pathways,
                                       enrichment = NULL,
                                       pathway_col = "pathway",
                                       fdr_col = "fdr",
                                       direction_col = "dir",
                                       weight_col = "weights") {
  out <- data.frame(pathway = pathways, stringsAsFactors = FALSE)
  out$fdr <- NA_real_
  out$dir <- NA_real_
  out$weights <- NA_real_

  if (is.null(enrichment)) {
    return(out)
  }
  if (!pathway_col %in% colnames(enrichment)) {
    stop("Missing enrichment pathway column: ", pathway_col)
  }

  idx <- match(pathways, enrichment[[pathway_col]])
  if (fdr_col %in% colnames(enrichment)) {
    out$fdr <- as.numeric(enrichment[[fdr_col]][idx])
  }
  if (direction_col %in% colnames(enrichment)) {
    out$dir <- as.numeric(enrichment[[direction_col]][idx])
  }
  if (weight_col %in% colnames(enrichment)) {
    out$weights <- as.numeric(enrichment[[weight_col]][idx])
  } else if (all(c(fdr_col, direction_col) %in% colnames(enrichment))) {
    out$weights <- -log10(pmax(out$fdr, .Machine$double.xmin)) * sign(out$dir)
  }

  out
}

write_json_like_config <- function(config, file) {
  lines <- json_value(config, indent = 0)
  writeLines(lines, file)
}

json_value <- function(x, indent = 0) {
  sp <- paste(rep(" ", indent), collapse = "")
  sp2 <- paste(rep(" ", indent + 2), collapse = "")
  if (is.list(x) && !is.null(names(x))) {
    fields <- vapply(names(x), function(name) {
      value <- json_value(x[[name]], indent = indent + 2)
      paste0(sp2, "\"", escape_json_string(name), "\": ", paste(value, collapse = "\n"))
    }, character(1))
    return(c("{", paste(fields, collapse = ",\n"), paste0(sp, "}")))
  }
  if (is.list(x)) {
    items <- vapply(x, function(item) paste(json_value(item, indent = indent + 2), collapse = "\n"), character(1))
    return(c("[", paste(paste0(sp2, items), collapse = ",\n"), paste0(sp, "]")))
  }
  if (length(x) != 1) {
    stop("JSON config values must be scalar or lists.")
  }
  if (is.numeric(x) || is.logical(x)) {
    return(tolower(as.character(x)))
  }
  paste0("\"", escape_json_string(as.character(x)), "\"")
}

escape_json_string <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  gsub("\"", "\\\\\"", x)
}

summarize_sumer_module_table <- function(module_table) {
  if (nrow(module_table) == 0) {
    return(data.frame())
  }
  modules <- sort(unique(module_table$module))
  do.call(rbind, lapply(modules, function(module_id) {
    x <- module_table[module_table$module == module_id, , drop = FALSE]
    weights <- if ("weights" %in% colnames(x)) x$weights else rep(NA_real_, nrow(x))
    fdr <- if ("fdr" %in% colnames(x)) x$fdr else rep(NA_real_, nrow(x))
    top_idx <- if (all(is.na(weights))) 1 else which.max(abs(weights))
    n_genes <- if ("n_genes" %in% colnames(x)) x$n_genes else rep(NA_integer_, nrow(x))
    data.frame(
      module = module_id,
      n_pathways = nrow(x),
      summed_pathway_genes = suppressWarnings(sum(n_genes, na.rm = TRUE)),
      top_pathway = x$pathway[top_idx],
      max_abs_weight = suppressWarnings(max(abs(weights), na.rm = TRUE)),
      min_fdr = suppressWarnings(min(fdr, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
}
