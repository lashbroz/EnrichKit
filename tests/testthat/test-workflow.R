test_that("GMT files can be read into pathway databases", {
  gmt <- file.path(tempdir(), "toy.gmt")
  writeLines(
    c(
      "MAPK\tna\tBRAF\tMAP2K1\tMAPK1",
      "PI3K\tna\tAKT1\tMTOR"
    ),
    gmt
  )

  sets <- read_gmt(gmt, universe = c("BRAF", "MAP2K1", "AKT1", "MTOR"), min_size = 2)
  expect_equal(names(sets), c("MAPK", "PI3K"))
  expect_equal(sets$MAPK, c("BRAF", "MAP2K1"))

  db <- make_pathway_db_from_gmt(gmt, database = "toy", min_size = 2)
  expect_s3_class(db, "EnrichKit_pathway_db")
  expect_equal(pathway_metadata(db)$database, c("toy", "toy"))
})

test_that("pathway background matching is auditable and orderable", {
  sets <- list(
    z_path = c("A", "B", "C", "Z"),
    a_path = c("A", "B"),
    dropped = c("Z")
  )
  db <- make_pathway_db(sets, database = "toy", min_size = 1)

  matched <- match_pathway_background(
    db,
    background = c("A", "B", "C"),
    min_size = 2,
    order_by = "pathway",
    warn = FALSE
  )

  expect_equal(names(as_gene_sets(matched)), c("a_path", "z_path"))
  summary <- pathway_matching_summary(matched)
  expect_true(all(c("original_n_genes", "matched_n_genes", "retained", "drop_reason") %in% colnames(summary)))
  expect_equal(summary$matched_n_genes[summary$pathway == "z_path"], 3)
  expect_false(summary$retained[summary$pathway == "dropped"])
  expect_equal(summary$drop_reason[summary$pathway == "dropped"], "below_min_size")
})

test_that("MSigDB collection metadata is available", {
  out <- msigdb_collections("human")
  expect_true(all(c("collection", "label", "gmt_file") %in% colnames(out)))
  expect_true("H" %in% out$collection)
  expect_true("C9" %in% out$collection)

  files <- msigdb_expected_files(c("H", "C2"), species = "human")
  expect_equal(names(files), c("H", "C2"))
  expect_true(grepl("v2026.1.Hs", files[["H"]]))
})

test_that("selected MSigDB GMT files can be combined", {
  h_file <- file.path(tempdir(), "h.all.v2026.1.Hs.symbols.gmt")
  c2_file <- file.path(tempdir(), "c2.all.v2026.1.Hs.symbols.gmt")
  writeLines("HALLMARK_MAPK\tna\tBRAF\tMAP2K1\tMAPK1", h_file)
  writeLines("KEGG_MAPK_SIGNALING_PATHWAY\tna\tBRAF\tMAP2K1\tMAPK3", c2_file)

  db <- build_msigdb_pathway_db(
    files = c(H = h_file, C2 = c2_file),
    species = "human",
    min_size = 2
  )

  meta <- pathway_metadata(db)
  expect_equal(nrow(meta), 2)
  expect_true(all(c("MSigDB_H_v2026.1.Hs", "MSigDB_C2_v2026.1.Hs") %in% meta$database))
})

test_that("KidsFirst default pathway database removes KEGG_MEDICUS and adds canonical KEGG", {
  hope_file <- file.path(tempdir(), "hope.gmt")
  kegg_file <- file.path(tempdir(), "kegg.gmt")
  writeLines(
    c(
      "KEGG_MEDICUS_BAD\tna\tA\tB\tC",
      "HOPE_PATH\tna\tA\tB\tC\tD\tE\tF"
    ),
    hope_file
  )
  writeLines("KEGG_MAPK_SIGNALING_PATHWAY\tna\tA\tB\tC\tD\tE\tF", kegg_file)

  db <- build_kfirst_default_pathway_db(
    hope_gmt = hope_file,
    canonical_kegg_gmt = kegg_file,
    universe = LETTERS[1:10],
    min_size = 2,
    max_size = 249
  )

  sets <- as_gene_sets(db)
  expect_false("KEGG_MEDICUS_BAD" %in% names(sets))
  expect_true("KEGG_MAPK_SIGNALING_PATHWAY" %in% names(sets))
})

test_that("redundancy reduction marks similar lower-ranked pathways", {
  sets <- list(
    A = c("G1", "G2", "G3"),
    B = c("G1", "G2", "G3", "G4"),
    C = c("X1", "X2", "X3")
  )
  results <- data.frame(pathway = c("A", "B", "C"), p = c(0.001, 0.01, 0.02))

  sim <- pathway_similarity(sets)
  expect_true(all(c("jaccard", "gene_separation") %in% colnames(sim)))

  reduced <- reduce_redundant_pathways(results, sets, jaccard_cutoff = 0.7)
  expect_false(reduced$redundant[reduced$pathway == "A"])
  expect_true(reduced$redundant[reduced$pathway == "B"])
  expect_equal(reduced$representative_pathway[reduced$pathway == "B"], "A")
})

test_that("cascade membership data tracks pathway members", {
  sets <- list(A = c("G1", "G2"), B = c("G2", "G3"))
  out <- cascade_membership_data(sets, pathways = c("A", "B"), members = c("G1", "G2", "G3"))

  expect_equal(nrow(out), 6)
  expect_true(out$present[out$pathway == "A" & out$member == "G1"])
  expect_false(out$present[out$pathway == "A" & out$member == "G3"])
})

test_that("cascade threshold consolidation retains pathways with new members", {
  sets <- list(
    P1 = c("A", "B", "C", "D"),
    P2 = c("A", "B", "C"),
    P3 = c("E", "F"),
    P4 = c("A", "E")
  )
  results <- data.frame(pathway = c("P1", "P2", "P3", "P4"), p = c(0.001, 0.01, 0.02, 0.03))

  out <- cascade_threshold_consolidation(
    results,
    sets,
    min_new_members = 2,
    row_order = "input"
  )

  expect_true("P1" %in% out$retained_pathways)
  expect_true("P3" %in% out$retained_pathways)
  expect_false(out$results$cascade_retained[out$results$pathway == "P2"])
  expect_equal(out$results$cascade_n_new_members[out$results$pathway == "P3"], 2)
})

test_that("supplement table has stable publication columns", {
  results <- data.frame(
    pathway = c("A", "B"),
    database = c("db", "db"),
    method = "wilcoxon",
    alternative = "greater",
    dir = c(1, -1),
    effect = c(2, -2),
    p = c(0.01, 0.2),
    fdr = c(0.02, 0.2),
    fdr_scope = "database",
    n_background = 100,
    n_pathway = c(5, 6),
    pathway_genes = c("G1;G2", "G3;G4")
  )

  out <- format_enrichment_supplement(results, analysis_name = "toy")
  expect_true(all(c("analysis", "pathway", "database", "p_value", "fdr", "genes") %in% colnames(out)))
  expect_equal(out$analysis, c("toy", "toy"))
  expect_equal(out$genes[1], "G1;G2")
})
