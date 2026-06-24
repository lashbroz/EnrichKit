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
