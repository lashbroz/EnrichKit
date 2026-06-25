test_that("prepare_sumer_input writes SUMER input files", {
  gene_sets <- list(
    PATH_A = c("A", "B", "C"),
    PATH_B = c("C", "D"),
    PATH_C = c("E", "F")
  )
  enrichment <- data.frame(
    pathway = c("PATH_A", "PATH_B", "MISSING"),
    fdr = c(0.01, 0.2, 0.03),
    dir = c(1, -1, 1)
  )
  prefix <- file.path(tempdir(), "toy_sumer")
  out <- prepare_sumer_input(enrichment, gene_sets, prefix)

  expect_true(file.exists(out$gmt_file))
  expect_true(file.exists(out$data_file))
  expect_equal(out$sumer_data$pathway, c("PATH_A", "PATH_B"))
  expect_gt(out$sumer_data$weights[1], 0)
  expect_lt(out$sumer_data$weights[2], 0)
})

test_that("prepare_sumer_input can select pathways for SUMER", {
  gene_sets <- list(
    PATH_A = c("A", "B", "C"),
    PATH_B = c("C", "D"),
    PATH_C = c("E", "F"),
    PATH_D = c("G", "H")
  )
  enrichment <- data.frame(
    pathway = c("PATH_A", "PATH_B", "PATH_C", "PATH_D"),
    fdr = c(0.01, 0.20, 0.03, 0.04),
    dir = c(1, -1, -1, 1)
  )
  prefix <- file.path(tempdir(), "toy_sumer_selected")
  summary_file <- paste0(prefix, "_selection_summary.tsv")

  out <- prepare_sumer_input(
    enrichment,
    gene_sets,
    prefix,
    fdr_threshold = 0.05,
    top_n = 2
  )

  expect_equal(out$sumer_data$pathway, c("PATH_A", "PATH_C"))
  expect_true(file.exists(summary_file))
  expect_equal(out$selection_summary_file, summary_file)
  summary <- utils::read.table(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expect_true("fdr_threshold_0.05" %in% summary$step)
  expect_true("top_n_2" %in% summary$step)
  expect_equal(summary$n_pathways[summary$step == "fdr_threshold_0.05"], 3)
  expect_equal(summary$n_pathways[summary$step == "top_n_2"], 2)

  score_file <- utils::read.table(out$data_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  expect_equal(ncol(score_file), 2)
})

test_that("get_sumer.data compatibility wrappers write SUMER input files", {
  gene_sets <- list(PATH_A = c("A", "B", "C"), PATH_B = c("C", "D"))
  enrichment <- data.frame(
    pathway = c("PATH_A", "PATH_B"),
    fdr = c(0.01, 0.2),
    dir = c(1, -1)
  )

  old_prefix <- file.path(tempdir(), "old_get_sumer")
  new_prefix <- file.path(tempdir(), "new_get_sumer")
  old <- get_sumer.data(enrichment, gene_sets, old_prefix)
  new <- get_sumer_data(enrichment, gene_sets, new_prefix)

  expect_true(file.exists(old$data_file))
  expect_true(file.exists(new$data_file))
  expect_equal(old$sumer_data$pathway, new$sumer_data$pathway)
  expect_equal(old$sumer_data$weights, new$sumer_data$weights)
})

test_that("SUMER config template and workflow prepare files", {
  gene_sets <- list(PATH_A = c("A", "B", "C"), PATH_B = c("C", "D"))
  enrichment <- data.frame(pathway = c("PATH_A", "PATH_B"), fdr = c(0.01, 0.2), dir = c(1, -1))
  prefix <- file.path(tempdir(), "workflow_sumer")
  config_file <- file.path(tempdir(), "workflow_sumer_config.json")

  job <- sumer_workflow(
    enrichment = enrichment,
    gene_sets = gene_sets,
    out_prefix = prefix,
    config_file = config_file,
    run = FALSE
  )

  expect_true(file.exists(job$prep$gmt_file))
  expect_true(file.exists(job$prep$data_file))
  expect_true(file.exists(job$config_file))
  config_text <- paste(readLines(job$config_file), collapse = "\n")
  expect_true(grepl("project", config_text))
  expect_true(grepl("top_num", config_text))
  expect_true(grepl("gmt_file", config_text))
  expect_true(grepl("score_file", config_text))
})

test_that("read_sumer_modules handles edgeless SUMER output", {
  node_file <- file.path(tempdir(), "ap_sumer_nodelist.txt")
  edge_file <- file.path(tempdir(), "ap_sumer_edgelist.txt")
  write.table(
    data.frame(name = c("score_sumer_assoc_PATH_A", "score_sumer_assoc_PATH_B")),
    file = node_file,
    quote = FALSE,
    row.names = FALSE
  )
  file.create(edge_file)

  gene_sets <- list(PATH_A = c("A", "B"), PATH_B = c("C", "D"))
  enrichment <- data.frame(pathway = c("PATH_A", "PATH_B"), fdr = c(0.01, 0.2), dir = c(1, -1))

  out <- read_sumer_modules(edge_file, node_file, gene_sets, enrichment)
  expect_equal(nrow(out$module_table), 2)
  expect_equal(nrow(out$module_summary), 2)
  expect_true(all(c("module", "pathway", "fdr", "weights") %in% colnames(out$module_table)))
})
