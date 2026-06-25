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
  collections, or the packaged KidsFirst `gosets.all` default.
- Preserve database/source labels and match every pathway to the interrogated
  genes for a specific analysis, with an audit trail of retained, trimmed, and
  dropped sets for cross-referencing analyses performed on different platforms
  or feature universes.

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
- Document how legacy HOPE/KidsFirst utilities map into EnrichKit functions.

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

## Core Workflow

Start with a named pathway list: names are pathway identifiers, values are gene
symbols. Small examples should still look like real pathway biology.

```r
library(EnrichKit)

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
pathway_db <- match_pathway_background(pathway_db, interrogated_genes)
```

Always use the analysis-specific background: interrogated proteins for protein
analyses, interrogated phosphosites for phosphosite analyses, and so on.

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
group by database label. This is useful when a small interrogated gene universe
causes many pathways to collapse onto the same few retained genes.

## MSigDB and KidsFirst Defaults

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

## Kids First Studywide Pathway Database And Generation

EnrichKit ships the current Kids First studywide `gosets.all` database as
package data. This is the database used for Kids First pathway enrichment unless
an analysis explicitly supplies a different pathway database.

```r
data("kfirst_gosets_all")
data("kfirst_gosets_source")
data("kfirst_gosets_metadata")
data("kfirst_gene_universe")

length(kfirst_gosets_all)
table(kfirst_gosets_source$source)
kfirst_gosets_metadata
```

The packaged object currently contains:

- 8,969 total retained pathways.
- 12,339 genes in the Kids First/interrogated gene universe.
- 8,785 pathways from `HOPE_pathway_database_without_KEGG_MEDICUS`.
- 184 pathways from `MSigDB_c2_cp_kegg_v7_canonical`.
- Inclusive pathway-size filtering of 6 to 249 matched genes for this packaged
  object.

For most analyses, use the accessor:

```r
gosets <- get_kfirst_gosets()

pathway_db <- get_kfirst_gosets(as_pathway_db = TRUE)
```

To construct the Kids First studywide database from scratch, start from the
explicit source GMT files and the measured/interrogated gene universe:

```r
kfirst_db <- build_kfirst_default_pathway_db(
  hope_gmt = "path/to/pathway_database_HOPE.gmt",
  canonical_kegg_gmt = "path/to/c2.cp.kegg.v7.0.symbols.gmt",
  universe = interrogated_genes,
  min_size = 6,
  max_size = 249
)
```

The current Kids First studywide database is built from these explicitly chosen
source databases:

- `HOPE_pathway_database_without_KEGG_MEDICUS`: HOPE pathway database after
  removing `KEGG_MEDICUS*` pathways.
- `MSigDB_c2_cp_kegg_v7_canonical`: canonical KEGG pathways from MSigDB C2
  canonical pathways, added so standard names such as
  `KEGG_MAPK_SIGNALING_PATHWAY` are available.

The retained HOPE-side pathway families are:

- `HALLMARK`
- `GOBP`
- `GOMF`
- `REACTOME`
- `BIOCARTA`
- `MITO3`

The HOPE-side `KEGG_MEDICUS*` pathways are deliberately excluded. Canonical KEGG
is instead supplied from the explicit MSigDB C2 CP KEGG GMT
(`c2.cp.kegg.v7.0.symbols.gmt`) so pathway names follow the standard MSigDB KEGG
convention.

The packaged data were generated from:

```r
kfirst_gosets_metadata$raw_sources
# HOPE_pathway_database: "pathway_database_HOPE.gmt"
# MSigDB_c2_cp_kegg_v7_canonical: "c2.cp.kegg.v7.0.symbols.gmt"
```

The generation script is included in the package repository:

```sh
Rscript data-raw/create_kfirst_gosets_all.R \
  path/to/pathway_database_HOPE.gmt \
  path/to/c2.cp.kegg.v7.0.symbols.gmt \
  path/to/kids_first_gene_universe.txt \
  data
```

The script applies the following steps:

```r
hope_sets <- read_gmt(hope_gmt)
kegg_sets <- read_gmt(canonical_kegg_gmt)
kfirst_gene_universe <- sort(unique(readLines(gene_universe_file)))

hope_no_medicus <- hope_sets[!grepl("^KEGG_MEDICUS", names(hope_sets))]
canonical_to_add <- kegg_sets[!names(kegg_sets) %in% names(hope_no_medicus)]
combined <- c(hope_no_medicus, canonical_to_add)

kfirst_gosets_all <- filter_gosets(
  combined,
  gene_universe = kfirst_gene_universe,
  min_genes = 6,
  max_genes = 249
)
```

`filter_gosets()` intersects each pathway with the Kids First gene universe and
retains pathways with inclusive matched size
`6 <= matched_n_genes <= 249`. The companion `kfirst_gosets_source` table stores
the final source label and matched gene count for every retained pathway.

The KidsFirst default can also load an external already-built GMT when a frozen
database artifact is the desired input:

```r
kfirst_gosets_provenance()

kfirst_pathway_database_components(
  "path/to/gosets_all_kfirst_source.tsv"
)

kfirst_db <- load_kfirst_gosets_gmt(
  file = "path/to/gosets_all_kfirst.gmt",
  source_table = "path/to/gosets_all_kfirst_source.tsv",
  universe = interrogated_genes
)
```

For methods sections or audit trails, `kfirst_gosets_provenance()` returns a
citation-ready provenance statement, the component table, and the size filters:

```r
prov <- kfirst_gosets_provenance()

cat(prov$text)
prov$components
prov$filters
```

The packaged data objects can also be reproduced with
`data-raw/create_kfirst_gosets_all.R`, which expects the HOPE pathway GMT, the
canonical MSigDB C2 KEGG GMT, a one-column Kids First gene-universe text file,
and an output directory.

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
  top_num = 100,
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
  top_num = 100,
  platform_name = "pathway_enrichment_top100",
  platform_abbr = "pathway",
  run = FALSE
)
```

You can combine both arguments. In that case EnrichKit first applies
`fdr_threshold`, orders the remaining pathways by `abs(weights)`, deduplicates
pathways if requested, and then applies `top_n`.

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

The workflow mirrors the old HOPE/KidsFirst pattern:

1. Write the pathway GMT and pathway-weight score file.
2. Write a `myconfig_*.json` style SUMER config.
3. Inspect/edit that config if needed.
4. Run SUMER manually from R.
5. Read the resulting modules back into R.

The default EnrichKit workflow stops before running SUMER:

```r
gene_sets <- as_gene_sets(pathway_db)

sumer_job <- sumer_workflow(
  enrichment = wilcox_res,
  gene_sets = gene_sets,
  out_prefix = "pathway_enrichment",
  weight_col = "signed.fdr",
  top_num = 100,
  platform_name = "pathway_enrichment",
  platform_abbr = "pathway",
  run = FALSE
)

sumer_job$prep$gmt_file
sumer_job$prep$data_file
sumer_job$config_file
```

This creates the same style of files used in the existing workflow:

- `pathway_enrichment_pathways.gmt`: pathway definitions used by SUMER.
- `pathway_enrichment_data.txt`: two-column pathway/weight file.
- `pathway_enrichment_sumer_config.json`: editable SUMER config.

For older HOPE/KidsFirst scripts, `get_sumer.data()` is provided as a
compatibility wrapper around `prepare_sumer_input()`. New code can use the more
explicit `prepare_sumer_input()` name or the snake-case alias
`get_sumer_data()`.

The config is intentionally written as a plain editable file. For multi-platform
SUMER analyses, add additional entries to the `data` array, each with its own
platform label, GMT file, and score file. See the SUMER GitHub repository for
the full expected config structure and options:
https://github.com/bzhanglab/sumer.

The config will look like:

```json
{
  "project": "pathway_enrichment_sumer_output",
  "top_num": 100,
  "data": [
    {
      "platform_name": "pathway_enrichment",
      "platform_abbr": "pathway",
      "gmt_file": "pathway_enrichment_pathways.gmt",
      "score_file": "pathway_enrichment_data.txt"
    }
  ]
}
```

Meaning:

- `project`: SUMER project/output label.
- `top_num`: number of top weighted pathways SUMER should use.
- `platform_name`: readable name for the analysis.
- `platform_abbr`: short prefix used in SUMER node labels.
- `gmt_file`: pathway membership file.
- `score_file`: pathway score/weight file.

Then run SUMER manually from R, matching the old workflow:

```r
library(sumer)
sumer("pathway_enrichment_sumer_config.json", "pathway_enrichment_sumer_output")
```

After SUMER finishes, read the module outputs:

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

If desired, EnrichKit can call SUMER for you with `run = TRUE`, but the manual
`sumer(config, output_name)` step is usually clearer because SUMER config details
can vary across installations.

For a fully scripted run:

```r
library(sumer)

sumer_job <- sumer_workflow(
  enrichment = wilcox_res,
  gene_sets = gene_sets,
  out_prefix = "pathway_enrichment",
  weight_col = "signed.fdr",
  top_num = 100,
  platform_name = "pathway_enrichment",
  platform_abbr = "pathway",
  run = TRUE,
  output_name = "pathway_enrichment_sumer_output",
  overwrite = TRUE
)

sumer_job$modules$module_table
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

- deeper HOPE_AYA script-by-script migration
- harmonized gene-set support, pending Weiping's relevant harmonization files
- polished pkgdown-style reference site
- broader plot theming for publication figures

## Testing

```r
pkgload::load_all(".")
testthat::test_dir("tests/testthat")
```
