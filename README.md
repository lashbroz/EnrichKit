# EnrichKit

`EnrichKit` is an R package for standardized pathway and enrichment analysis
across CPTAC/Kids First-style proteogenomic projects. The goal is to use one explicit,
auditable system for pathway database construction, background matching,
enrichment testing, pathway consolidation, SUMER interoperability, and
publication-ready supplementary tables.

The package is designed for reproducible, interoperable cross-study analysis:
explicit inputs, explicit backgrounds, stable result columns, database-aware FDR
correction, and clear database provenance.

## What EnrichKit Does

**Prepare pathway databases**

- Build gene-set databases from named lists, GMT files, selected MSigDB
  collections, or the packaged Kids First pathway database.
- Preserve database/source labels and match every pathway to the measured assay
  background for a specific analysis, with an audit trail of retained, trimmed,
  and dropped sets for cross-referencing analyses performed on different
  platforms or feature sets.

**Run enrichment analyses**

- Test hit lists with Fisher/exact tests and ranked or scored features with
  Wilcoxon-based methods.
- Support one-sided, two-sided, score-thresholded, and rank-shift Wilcoxon
  analyses.
- Apply either database-wise or global FDR correction and return signed
  statistics for ranking, heatmaps, and SUMER weights.

**Consolidate pathway results**

- Prepare SUMER GMT, score, and config files; run SUMER when available; and read
  SUMER modules back into R.
- Provide transparent redundancy filters and cascade visualizations for tracking
  which genes/features are retained or lost during pathway consolidation.
- Keep original pathway names after consolidation through representative-to-all
  pathway crosswalks, so equivalent pathway signals can still be compared across
  analyses that used different underlying gene sets.

**Generate manuscript-ready outputs**

- Produce stable enrichment-result and supplementary-table formats with pathway
  labels, database labels, counts, p-values, FDR values, effect direction,
  method metadata, and background size.
- Document how legacy project utilities map into EnrichKit functions.

## Installation

From GitHub:

```r
install.packages("devtools")
devtools::install_github("lashbroz/EnrichKit")
```

From a local checkout:

```r
devtools::install(".")
```

During development, load directly:

```r
pkgload::load_all(".")
```

Optional SUMER support uses the Zhang lab SUMER package. Install it directly
from GitHub:

```r
devtools::install_github("bzhanglab/sumer")
library(sumer)

exists("sumer", mode = "function")
# TRUE
```

## Pathway Database Generation

```r
library(EnrichKit)
```

The examples below use `interrogated_genes` for the measured assay background:
the genes or features that were actually measured and eligible for testing in
the analysis. For a protein matrix, this is usually the protein/gene identifiers
represented in the matrix rows.

### MSigDB and KidsFirst Defaults

Use `msigdb_collections()` to see the supported current MSigDB collection menu
and expected GMT file names:

```r
msigdb_collections("human")
msigdb_expected_files(c("H", "C2", "C5"), species = "human")
```

After downloading GMT files from MSigDB, combine the selected collections:

```r
msigdb_db <- build_msigdb_pathway_db(
  files = c(
    H = "~/Downloads/h.all.v2026.1.Hs.symbols.gmt",
    C2 = "~/Downloads/c2.all.v2026.1.Hs.symbols.gmt",
    C5 = "~/Downloads/c5.all.v2026.1.Hs.symbols.gmt"
  ),
  species = "human",
  universe = interrogated_genes,
  min_size = 5,
  max_size = 250
)
```

### Kids First Studywide Pathway Database And Generation

EnrichKit ships the current Kids First studywide pathway database as package
data. This is the database used for Kids First pathway enrichment unless an
analysis explicitly supplies a different pathway database.

```r
data("kfirst_gosets_all")
data("kfirst_gene_universe")

length(kfirst_gosets_all)
# Measured assay background used to build the packaged database.
length(kfirst_gene_universe)
```

The packaged object currently contains:

- 9,776 total retained pathways.
- 12,339 genes in the Kids First measured assay background.
- Inclusive pathway-size filtering of 5 to 250 matched genes for this packaged
  object.

For most analyses, use the accessor:

```r
pathway_db <- get_kfirst_gosets(as_pathway_db = TRUE)
```

To construct a Kids First studywide database from scratch, start from explicit
MSigDB GMT files and the measured assay background:

```r
msigdb_db <- build_msigdb_pathway_db(
  files = c(
    H = "path/to/h.all.v2026.1.Hs.symbols.gmt",
    C2 = "path/to/c2.all.v2026.1.Hs.symbols.gmt",
    C5 = "path/to/c5.all.v2026.1.Hs.symbols.gmt"
  ),
  species = "human",
  universe = interrogated_genes,
  min_size = 5,
  max_size = 250
)
```

The recommended MSigDB sources for broad proteogenomic pathway screening are:

- `HALLMARK`
- `C2:CP:REACTOME`
- `C2:CP:BIOCARTA`
- `C2:CP:KEGG`
- `C5:GO:BP`
- `C5:GO:MF`

The same filtering logic applies to any explicit GMT source:

```r
pathway_db <- match_pathway_background(
  msigdb_db,
  interrogated_genes,
  min_size = 5,
  max_size = 250
)
```

Background matching intersects each pathway with the measured assay background
and retains pathways with inclusive matched size
`5 <= matched_n_genes <= 250`.

The KidsFirst default can also load an external already-built GMT when a frozen
database artifact is the desired input:

```r
kfirst_db <- load_kfirst_gosets_gmt(
  file = "path/to/gosets_all_kfirst.gmt",
  source_table = "path/to/gosets_all_kfirst_source.tsv",
  universe = interrogated_genes
)
```

For reproducible projects, keep the exact GMT files, MSigDB version labels, and
the one-column measured assay background file with the analysis outputs.

## Core Enrichment Workflow

Once a pathway database has been built or loaded, match it to the
analysis-specific background before running enrichment: interrogated proteins
for protein analyses, interrogated phosphosites for phosphosite analyses, and so
on.

If starting from a small named list, names are pathway identifiers and values are
gene symbols. Small examples should still look like real pathway biology.

```r
gene_sets <- list(
  KEGG_MAPK_SIGNALING_PATHWAY = c("BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RAF1", "KRAS", "NRAS"),
  BIOCARTA_MTOR_PATHWAY = c("PIK3CA", "PIK3R1", "AKT1", "AKT2", "MTOR", "RPTOR", "RICTOR"),
  REACTOME_CELL_CYCLE_CHECKPOINTS = c("CDK1", "CDK2", "RB1", "CHEK1", "CHEK2", "CCNB1", "CCNE1")
  # ...
)

pathway_db <- make_pathway_db(
  gene_sets,
  database = c(
    KEGG_MAPK_SIGNALING_PATHWAY = "KEGG",
    BIOCARTA_MTOR_PATHWAY = "BIOCARTA",
    REACTOME_CELL_CYCLE_CHECKPOINTS = "REACTOME"
  ),
  source = "example_pathway_database",
  min_size = 5,
  max_size = 250
)

interrogated_genes <- c(
  "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RAF1", "KRAS", "NRAS",
  "PIK3CA", "PIK3R1", "AKT1", "AKT2", "MTOR", "RPTOR", "RICTOR",
  "CDK1", "CDK2", "RB1", "CHEK1", "CHEK2", "CCNB1", "CCNE1", "TP53"
)
```

`match_pathway_background()` is intentionally explicit because it changes the
effective pathway database. It intersects each pathway with the interrogated
background, drops pathways whose matched gene count falls outside the inclusive
`min_size <= matched_n_genes <= max_size` range, sorts members within each
pathway, and records what happened:

```r
pathway_db <- match_pathway_background(
  pathway_db,
  interrogated_genes,
  min_size = 5,
  max_size = 250,
  order_by = "input"
)

pathway_matching_summary(pathway_db)
```

The default `order_by = "input"` preserves the source database order. Use
`order_by = "pathway"` for alphabetical order, or `order_by = "database"` to
group by database label. This is useful when a small measured assay background
causes many pathways to collapse onto the same few retained genes.

## Fisher Enrichment

```r
hits <- c("BRAF", "MAP2K1", "MAPK1")

fisher_res <- fisher_enrichment(
  hits = hits,
  pathway_db = pathway_db,
  background = interrogated_genes,
  alternative = "greater"
)
```

`alternative` is explicit:

- `"greater"` tests pathway over-representation among hits.
- `"less"` tests pathway under-representation among hits.
- `"two.sided"` tests either direction.

By default, FDR is adjusted within each pathway database:

```r
unique(fisher_res$fdr_scope)
# "database"
```

Use `fdr_by_database = FALSE` for one global adjustment across all pathways.

## Wilcoxon Enrichment

```r
scores <- c(
  BRAF = 2.4,
  MAP2K1 = 1.7,
  MAPK1 = 1.3,
  MAPK3 = 1.1,
  AKT1 = -0.3,
  MTOR = -0.8,
  TP53 = 0.2
)

wilcox_res <- wilcox_enrichment(
  feature_scores = scores,
  pathway_db = pathway_db,
  background = interrogated_genes,
  alternative = "greater"
)
```

For score-thresholded Wilcoxon enrichment:

```r
thresholded_res <- thresholded_wilcox_enrichment(
  feature_scores = scores,
  pathway_db = pathway_db,
  threshold = 1,
  alternative = "greater"
)
```

For the rank-shift Wilcoxon method from Xiaoyu Song:

```r
rank_shift_res <- rank_shift_wilcox_enrichment(
  feature_scores = scores,
  pathway_db = pathway_db,
  background = interrogated_genes,
  shift_fraction = 0.20,
  direction = "both"
)
```

This method converts scores to ranks and performs one-sided tests separately
for enriched and depleted pathways. For enrichment, non-pathway ranks are
shifted upward by `round(n_features * shift_fraction)` before testing whether
pathway ranks are still greater. For depletion, non-pathway ranks are shifted
downward before testing whether pathway ranks are still smaller.

Use this method when ordinary Wilcoxon enrichment is too permissive and the
goal is to identify pathways whose signal remains convincing after imposing a
prespecified rank margin against the background. It is best treated as a
robustness or stringency analysis: report the `shift_fraction`, and consider
checking whether conclusions are stable across nearby values such as `0.10`,
`0.20`, and `0.30`.

## Result Columns

Core enrichment outputs include:

- `pathway`
- `database`
- `n_background`
- pathway/hit counts
- `effect` or `odds_ratio`
- `dir`
- `p`
- `fdr`
- `signed.p`
- `signed.fdr`
- `alternative`
- `fdr_scope`
- `method`

The `signed.*` columns are signed `-log10(p)` values, using `dir` for the sign.

## Pathway Consolidation With SUMER

Enrichment often returns many overlapping pathways: related GO terms, Reactome
branches, KEGG/BioCarta variants, and near-duplicate signatures. SUMER is the
primary EnrichKit route for pathway consolidation when the goal is to summarize
many overlapping enrichment calls into network modules of related pathways.
SUMER was developed by the Zhang lab and remains an external tool:
https://github.com/bzhanglab/sumer. EnrichKit does not reimplement SUMER; it
handles the file preparation and round-trip workflow.

### Select pathways before SUMER

Do not assume that every tested pathway should be handed to SUMER. Large pathway
databases can return thousands of nominal or FDR-significant sets, especially
for broad proteomic contrasts. Passing all tested pathways to SUMER can make the
run slow and can obscure the biology with weak, redundant, or effectively
zero-weight nodes.

`prepare_sumer_input()` and `sumer_workflow()` can do this pathway selection
directly. Common choices are:

- pathways passing an analysis-defined threshold, such as
  `fdr_threshold = 0.10`;
- the top `N` pathways by absolute SUMER weight, such as `top_n = 100`;
- pathways from a more stringent method, such as
  `rank_shift_wilcox_enrichment()`;
- Fisher enrichment results from a high-confidence feature hit list, such as
  genes/proteins with feature-level `fdr < 0.10`, tested separately for up and
  down directions.

For example, to hand SUMER only FDR-significant pathways:

```r
sumer_job <- sumer_workflow(
  enrichment = wilcox_res,
  gene_sets = gene_sets,
  out_prefix = "pathway_enrichment_fdr10",
  weight_col = "signed.fdr",
  fdr_threshold = 0.10,
  platform_name = "pathway_enrichment_fdr10",
  platform_abbr = "pathway",
  run = FALSE
)

sumer_job$prep$selection_summary
sumer_job$prep$selection_summary_file
```

Or, keep only the strongest weighted pathways:

```r
sumer_job <- sumer_workflow(
  enrichment = wilcox_res,
  gene_sets = gene_sets,
  out_prefix = "pathway_enrichment_top100",
  weight_col = "signed.fdr",
  top_n = 100,
  platform_name = "pathway_enrichment_top100",
  platform_abbr = "pathway",
  run = FALSE
)
```

Treat `top_n` and `top_num` as different layers. `top_n` is an EnrichKit
pre-filter: it controls how many pathways are written into the SUMER input
files. `top_num` is a SUMER config field: it controls how SUMER uses the already
written pathway files. In most analyses, make the top-pathway decision once with
`top_n` or with an FDR threshold, then leave `top_num` at its default unless you
are intentionally editing the SUMER config.

You can combine `fdr_threshold` and `top_n` when that is the analysis design. In
that case EnrichKit first applies `fdr_threshold`, orders the remaining pathways
by `abs(weights)`, deduplicates pathways if requested, and then applies `top_n`.

Every SUMER input preparation writes a pathway-selection summary TSV by default:

```r
sumer_job$prep$selection_summary_file
# "pathway_enrichment_top100_selection_summary.tsv"
```

The summary records how many pathways remain after each step, including the FDR
threshold and/or top-N filter. This makes it easy to report exactly how many
pathways were handed to SUMER.

Treat this filtering as part of the analysis design and report it. SUMER is a
consolidation step, not a replacement for choosing a defensible enrichment
threshold or ranking rule.

Install SUMER before running SUMER jobs:

```r
devtools::install_github("bzhanglab/sumer")
library(sumer)

stopifnot(exists("sumer", mode = "function"))
```

The default EnrichKit workflow keeps SUMER file generation and configuration
under the hood:

1. Select pathways to pass to SUMER.
2. Write the pathway and weight files SUMER needs.
3. Write the SUMER configuration needed for the run.
4. Optionally run SUMER.
5. Read the resulting modules back into R.

For a fully scripted run:

```r
library(sumer)

gene_sets <- as_gene_sets(pathway_db)

sumer_job <- sumer_workflow(
  enrichment = wilcox_res,
  gene_sets = gene_sets,
  out_prefix = "pathway_enrichment",
  weight_col = "signed.fdr",
  platform_name = "pathway_enrichment",
  platform_abbr = "pathway",
  run = TRUE,
  output_name = "pathway_enrichment_sumer_output",
  overwrite = TRUE
)

sumer_job$modules$module_table
sumer_job$modules$module_summary
```

If you want to prepare files but run SUMER yourself, set `run = FALSE`:

```r
sumer_job <- sumer_workflow(
  enrichment = wilcox_res,
  gene_sets = gene_sets,
  out_prefix = "pathway_enrichment",
  weight_col = "signed.fdr",
  platform_name = "pathway_enrichment",
  platform_abbr = "pathway",
  run = FALSE
)

sumer_job$prep$selection_summary
sumer_job$prep$gmt_file
sumer_job$prep$data_file
```

EnrichKit writes the files SUMER needs. The configuration file is also written,
but most users can treat it as internal unless they are customizing SUMER:

- `pathway_enrichment_pathways.gmt`: pathway definitions used by SUMER.
- `pathway_enrichment_data.txt`: two-column pathway/weight file.

### Custom SUMER Configuration

Most single-analysis workflows should not need to edit the SUMER config. Treat
the JSON file as an implementation detail unless you are deliberately using
custom SUMER options.

SUMER can support multiple platforms, but that is not EnrichKit's native default
setup. For multi-platform SUMER analyses, provide a custom SUMER configuration
with one data entry per platform, each with its own platform label, GMT file, and
score file. See the SUMER repository for the full expected config structure:
https://github.com/bzhanglab/sumer.

If SUMER was run outside EnrichKit, read module outputs explicitly:

```r
modules <- read_sumer_modules(
  edge_file = "pathway_enrichment_sumer_output/ap_sumer_edgelist.txt",
  node_file = "pathway_enrichment_sumer_output/ap_sumer_nodelist.txt",
  gene_sets = gene_sets,
  enrichment = wilcox_res
)

modules$module_table
modules$module_summary
```

For cross-analysis comparison, keep a pathway-key table in addition to the
consolidated display table. For example, if a collaborator runs an analysis on a
subset of samples and you want to ask whether the same pathways appear in your
protein or RNA analysis, use the original pathway names as stable keys:

```r
protein_reduced <- reduce_redundant_pathways(protein_res, gene_sets)
subset_reduced <- reduce_redundant_pathways(collaborator_subset_res, gene_sets)

protein_map <- pathway_consolidation_map(protein_reduced)
subset_map <- pathway_consolidation_map(subset_reduced)

pathway_keys <- cross_reference_pathways(
  list(
    protein = protein_reduced,
    collaborator_subset = subset_reduced
  )
)

pathway_keys$wide
```

The `wide` table has one row per original pathway key and records whether that
pathway was present in each analysis, which representative pathway it was
consolidated under, and the corresponding FDR/direction when available. This is
the table to merge across platforms or append as a publication audit trail.

## Transparent Cascade Consolidation

For publication figures and ordered enrichment tables, EnrichKit also includes a
transparent cascade threshold method for tracking gene/feature content through
the consolidation process. Walking down a ranked pathway list, a pathway is
retained only if it introduces at least a specified number of new genes/features
beyond the pathways already retained. The output makes it clear which genes are
newly contributed, already covered, or lost when redundant pathways are removed.

```r
cascade <- cascade_threshold_consolidation(
  wilcox_res,
  gene_sets,
  min_new_members = 5,
  row_order = "input"
)

cascade$results
plot_cascade_threshold_consolidation(cascade)
```

Use SUMER when you want network-based module discovery. Use cascade threshold
consolidation when you need a simple, ordered, publication-facing pathway list
and an auditable record of gene-content retained or lost during consolidation.

## Development Status

Implemented:

- pathway database construction
- background matching
- pathway gene counts
- Fisher enrichment
- Wilcoxon enrichment
- thresholded Wilcoxon enrichment
- one-tailed and two-tailed test specificity
- database-wise or global FDR adjustment
- signed p-value helpers
- SUMER input preparation
- SUMER execution wrapper
- SUMER module extraction
- MSigDB/GMT import helpers
- redundancy reduction
- cascade/pathway membership data and base plot
- manuscript supplementary table formatter
- minimal vignette
- migration map in `data-raw/migration_map.tsv`
- unit tests for the core pieces

Still in progress:

- deeper legacy script-by-script migration
- harmonized gene-set support, pending Weiping's relevant harmonization files
- polished pkgdown-style reference site
- broader plot theming for publication figures

## Testing

```r
pkgload::load_all(".")
testthat::test_dir("tests/testthat")
```
