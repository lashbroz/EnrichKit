test_that("clean_feature_sets filters sets", {
  sets <- list(a = c("A", "A", "B"), b = "C")
  out <- clean_feature_sets(sets, universe = c("A", "B"), min_size = 2)
  expect_equal(names(out), "a")
  expect_equal(out$a, c("A", "B"))
})

test_that("calc_ora returns enrichment table", {
  sets <- list(mapk = c("MAPK1", "MAPK3", "BRAF"), pi3k = c("AKT1", "MTOR"))
  hits <- c("MAPK1", "BRAF")
  universe <- unique(c(unlist(sets), hits, "TP53"))

  out <- calc_ora(hits, sets, universe = universe)
  expect_true(all(c("pathway", "p", "fdr", "overlap_genes") %in% colnames(out)))
  expect_equal(out$n_hits[1], 2)
})

test_that("make_pathway_db preserves named database labels", {
  sets <- list(mapk = c("MAPK1", "MAPK3"), pi3k = c("AKT1", "MTOR"))
  db <- make_pathway_db(
    sets,
    database = c(mapk = "kinase", pi3k = "kinase"),
    source = "toy"
  )
  expect_equal(pathway_metadata(db)$database, c("kinase", "kinase"))
})

test_that("Fisher enrichment exposes one- and two-tailed alternatives", {
  sets <- list(path = c("A", "B"), other = c("E", "F"))
  hits <- c("A", "B")
  background <- LETTERS[1:6]

  greater <- fisher_enrichment(hits, sets, background = background, alternative = "greater")
  less <- fisher_enrichment(hits, sets, background = background, alternative = "less")
  two_sided <- fisher_enrichment(hits, sets, background = background, alternative = "two.sided")

  path_greater <- greater[greater$pathway == "path", ]
  path_less <- less[less$pathway == "path", ]
  path_two_sided <- two_sided[two_sided$pathway == "path", ]

  expect_equal(path_greater$alternative, "greater")
  expect_equal(path_less$alternative, "less")
  expect_equal(path_two_sided$alternative, "two.sided")
  expect_lt(path_greater$p, path_less$p)
  expect_gt(path_greater$dir, 0)
})

test_that("Wilcoxon enrichment exposes directional alternatives", {
  scores <- c(A = 5, B = 4, C = 0, D = 0, E = -1, F = -2)
  sets <- list(high_path = c("A", "B"), low_path = c("E", "F"))

  greater <- wilcox_enrichment(scores, sets, alternative = "greater")
  less <- wilcox_enrichment(scores, sets, alternative = "less")

  high_greater <- greater[greater$pathway == "high_path", ]
  high_less <- less[less$pathway == "high_path", ]
  low_less <- less[less$pathway == "low_path", ]

  expect_equal(high_greater$alternative, "greater")
  expect_lt(high_greater$p, high_less$p)
  expect_gt(high_greater$effect, 0)
  expect_lt(low_less$effect, 0)
})

test_that("thresholded Wilcoxon records threshold metadata", {
  scores <- c(A = 5, B = 4, C = 0.2, D = 0.1, E = -3, F = -2)
  sets <- list(high_path = c("A", "B"), low_path = c("E", "F"))

  out <- thresholded_wilcox_enrichment(
    scores,
    sets,
    threshold = 2,
    alternative = "two.sided"
  )

  expect_true(all(out$method == "thresholded_wilcoxon"))
  expect_true(all(out$threshold == 2))
  expect_true(all(out$n_threshold_background == 4))
})

test_that("rank-shift Wilcoxon implements Xiaoyu Song direction-specific tests", {
  scores <- c(A = 6, B = 5, C = 4, D = 3, E = 2, F = 1)
  sets <- list(high_path = c("A", "B"), low_path = c("E", "F"))

  out <- rank_shift_wilcox_enrichment(
    scores,
    sets,
    shift_fraction = 0.20,
    direction = "both"
  )

  expect_true(all(c("enriched", "depleted") %in% out$direction_test))
  expect_true(all(out$method == "rank_shift_wilcoxon_xiaoyu_song"))
  expect_true(all(out$rank_shift == round(length(scores) * 0.20)))
  expect_equal(
    out$alternative[out$pathway == "high_path" & out$direction_test == "enriched"],
    "greater"
  )
  expect_equal(
    out$alternative[out$pathway == "low_path" & out$direction_test == "depleted"],
    "less"
  )
})

test_that("enrichment defaults to database-wise FDR adjustment", {
  sets <- list(
    db1_a = c("A", "B"),
    db1_b = c("A", "C"),
    db2_a = c("A", "B", "D")
  )
  db <- make_pathway_db(
    sets,
    database = c(db1_a = "db1", db1_b = "db1", db2_a = "db2")
  )

  out_db <- fisher_enrichment(
    hits = c("A", "B"),
    pathway_db = db,
    background = LETTERS[1:6],
    alternative = "greater",
    fdr_by_database = TRUE
  )
  out_global <- fisher_enrichment(
    hits = c("A", "B"),
    pathway_db = db,
    background = LETTERS[1:6],
    alternative = "greater",
    fdr_by_database = FALSE
  )

  db2_fdr <- out_db$fdr[out_db$pathway == "db2_a"]
  global_db2_fdr <- out_global$fdr[out_global$pathway == "db2_a"]

  expect_equal(unique(out_db$fdr_scope), "database")
  expect_equal(unique(out_global$fdr_scope), "global")
  expect_equal(db2_fdr, out_db$p[out_db$pathway == "db2_a"])
  expect_gt(global_db2_fdr, db2_fdr)
})
