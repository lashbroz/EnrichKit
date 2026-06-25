#!/usr/bin/env Rscript

## Rebuild the packaged Kids First gosets.all object from source GMT files.
##
## Inputs:
##   1. HOPE pathway database GMT.
##   2. Canonical MSigDB C2 KEGG GMT.
##   3. Kids First/interrogated gene universe, one gene symbol per line.
##
## Provenance:
##   - Source 1: HOPE pathway database GMT.
##   - Source 2: canonical MSigDB C2 KEGG GMT.
##   - Remove KEGG_MEDICUS* pathways from the HOPE pathway database.
##   - Add canonical MSigDB C2 KEGG pathways not already present.
##   - Intersect each pathway with the Kids First gene universe.
##   - Retain pathways with inclusive matched size 5 to 250 genes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "Usage: Rscript data-raw/create_kfirst_gosets_all.R ",
    "<hope_gmt> <canonical_kegg_gmt> <gene_universe_txt> <out_dir>",
    call. = FALSE
  )
}

hope_gmt <- args[[1]]
canonical_kegg_gmt <- args[[2]]
gene_universe_file <- args[[3]]
out_dir <- args[[4]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

min_genes <- 5
max_genes <- 250

read_gmt <- function(file) {
  lines <- readLines(file, warn = FALSE)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  out <- lapply(parts, function(x) {
    genes <- unique(x[-c(1, 2)])
    genes[!is.na(genes) & nzchar(genes)]
  })
  names(out) <- vapply(parts, function(x) x[[1]], character(1))
  out
}

filter_gosets <- function(gosets, gene_universe, min_genes, max_genes) {
  out <- lapply(gosets, function(x) sort(unique(intersect(x, gene_universe))))
  n_genes <- lengths(out)
  out[n_genes >= min_genes & n_genes <= max_genes]
}

hope_sets <- read_gmt(hope_gmt)
kegg_sets <- read_gmt(canonical_kegg_gmt)
kfirst_gene_universe <- sort(unique(readLines(gene_universe_file, warn = FALSE)))
kfirst_gene_universe <- kfirst_gene_universe[nzchar(kfirst_gene_universe)]

hope_no_medicus <- hope_sets[!grepl("^KEGG_MEDICUS", names(hope_sets))]
canonical_to_add <- kegg_sets[!names(kegg_sets) %in% names(hope_no_medicus)]
combined <- c(hope_no_medicus, canonical_to_add)

source_lookup <- c(
  stats::setNames(
    rep("HOPE_pathway_database_without_KEGG_MEDICUS", length(hope_no_medicus)),
    names(hope_no_medicus)
  ),
  stats::setNames(
    rep("MSigDB_c2_cp_kegg_v7_canonical", length(canonical_to_add)),
    names(canonical_to_add)
  )
)

kfirst_gosets_all <- filter_gosets(
  combined,
  gene_universe = kfirst_gene_universe,
  min_genes = min_genes,
  max_genes = max_genes
)

kfirst_gosets_source <- data.frame(
  pathway = names(kfirst_gosets_all),
  source = unname(source_lookup[names(kfirst_gosets_all)]),
  n_genes = lengths(kfirst_gosets_all),
  stringsAsFactors = FALSE
)

kfirst_gosets_metadata <- list(
  name = "Kids First gosets.all",
  min_genes = min_genes,
  max_genes = max_genes,
  n_pathways = length(kfirst_gosets_all),
  n_gene_universe = length(kfirst_gene_universe),
  source_counts = as.list(table(kfirst_gosets_source$source)),
  raw_sources = c(
    HOPE_pathway_database = basename(hope_gmt),
    MSigDB_c2_cp_kegg_v7_canonical = basename(canonical_kegg_gmt)
  ),
  removed = "KEGG_MEDICUS* pathways removed from HOPE pathway database",
  added = "Canonical MSigDB C2 KEGG pathways not already present in HOPE after KEGG_MEDICUS removal"
)

save(kfirst_gosets_all, file = file.path(out_dir, "kfirst_gosets_all.rda"), compress = "xz")
save(kfirst_gosets_source, file = file.path(out_dir, "kfirst_gosets_source.rda"), compress = "xz")
save(kfirst_gene_universe, file = file.path(out_dir, "kfirst_gene_universe.rda"), compress = "xz")
save(kfirst_gosets_metadata, file = file.path(out_dir, "kfirst_gosets_metadata.rda"), compress = "xz")

utils::write.table(
  kfirst_gosets_source,
  file = file.path(out_dir, "kfirst_gosets_source.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Final filtered gene sets: ", length(kfirst_gosets_all))
message("Source counts:")
print(table(kfirst_gosets_source$source))
