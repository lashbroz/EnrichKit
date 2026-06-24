# EnrichKit

`EnrichKit` is an R package for standardized pathway and enrichment analysis
across CPTAC/Kids First-style proteogenomic projects. The goal is to use one explicit,
auditable system for pathway database construction, background matching,
enrichment testing, pathway consolidation, SUMER interoperability, and
publication-ready supplementary tables.

The package is designed for reproducible, interoperable cross-study analysis:
explicit inputs, explicit backgrounds, stable result columns, database-aware FDR
correction, and no hidden project paths.

## What EnrichKit Does

**1. Prepare standardized pathway databases for consistent use across a study**

- Builds pathway/gene-set databases from named lists, GMT files, selected MSigDB
  collections, or KidsFirst defaults.
- Reports supported MSigDB collection names and expected GMT filenames so users
  can choose databases explicitly.
- Combines selected databases while preserving database labels such as KEGG,
  BIOCARTA, REACTOME, GO, Hallmark, or KidsFirst source labels.
- Matches pathway databases to the interrogated gene/feature universe for a
  given analysis and records an audit table of retained, trimmed, and dropped
  pathways.
- Supports consistent database and background handling across analyses where
  pathway/gene-set interpretation is appropriate.
- Facilitates reproducibility and inter-study analytical interoperability by
  making database provenance, background matching, size filters, and pathway
  membership changes explicit.

**2. Implement pathway enrichment methods**

- Fisher/exact-test enrichment for hit lists, with explicit `greater`, `less`,
  and `two.sided` alternatives.
- Wilcoxon pathway enrichment for ranked or scored features.
- Score-thresholded Wilcoxon enrichment.
- Rank-shift Wilcoxon enrichment, a method from Xiaoyu Song.
- Database-wise or global FDR adjustment.
- Signed `-log10(p)` and signed FDR values for heatmaps, ranking, and SUMER
  weights.

**3. Interoperate with SUMER for pathway-module summaries**

- Writes pathway GMT files and pathway-weight files for SUMER.
- Writes editable `myconfig_*.json`-style SUMER config files.
- Supports the existing manual workflow: prepare files in R, inspect/reuse the
  config, run `sumer(config, output_name)`, and read modules back into R.
- Reads SUMER edge/node outputs into module tables and module summaries.

**4. Consolidate redundant enrichment results**

- Redundancy filters based on pathway overlap or fixed gene-separation rules.
- Transparent cascade consolidation for ordered pathway lists, with membership
  visualizations showing which genes/features are newly contributed by each
  retained pathway.

**5. Produce publication-oriented outputs**

- Stable enrichment-result columns for manuscripts and supplements.
- Supplementary table formatting with pathway labels, database labels, counts,
  p-values, FDR values, direction/effect, method, and background size.
- Migration-map documentation for how HOPE/KidsFirst-style utilities map into
  EnrichKit functions.

## Installation

From GitHub:

```r
install.packages("devtools")
devtools::install_github("lashbn01/EnrichKit")
```

From a local checkout:

```r
devtools::install(".")
```

During development, load directly:

```r
pkgload::load_all(".")
```

## Core Workflow

```r
library(EnrichKit)

gene_sets <- list(
  KEGG_MAPK_SIGNALING_PATHWAY = c("BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3"),
  BIOCARTA_MTOR_PATHWAY = c("PIK3CA", "AKT1", "MTOR"),
  REACTOME_CELL_CYCLE_CHECKPOINTS = c("CDK1", "CDK2", "RB1")
)

pathway_db <- make_pathway_db(
  gene_sets,
  database = c(
    KEGG_MAPK_SIGNALING_PATHWAY = "KEGG",
    BIOCARTA_MTOR_PATHWAY = "BIOCARTA",
    REACTOME_CELL_CYCLE_CHECKPOINTS = "REACTOME"
  ),
  source = "example_pathway_database",
  min_size = 2
)

interrogated_genes <- c("BRAF", "MAP2K1", "MAPK1", "MAPK3", "AKT1", "MTOR", "TP53", "CDK1")
pathway_db <- match_pathway_background(pathway_db, interrogated_genes)
```

Always use the analysis-specific background: interrogated proteins for protein
analyses, interrogated phosphosites for phosphosite analyses, and so on.

`match_pathway_background()` is intentionally explicit because it changes the
effective pathway database. It intersects each pathway with the interrogated
background, drops pathways outside the size limits, sorts members within each
pathway, and records what happened:

```r
pathway_db <- match_pathway_background(
  pathway_db,
  interrogated_genes,
  min_size = 6,
  max_size = 249,
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
  min_size = 6,
  max_size = 249
)
```

EnrichKit ships the current KidsFirst `gosets.all` database as package data:

```r
data("kfirst_gosets_all")
data("kfirst_gosets_source")

length(kfirst_gosets_all)
table(kfirst_gosets_source$source)
```

For most analyses, use the accessor:

```r
gosets <- get_kfirst_gosets()

pathway_db <- get_kfirst_gosets(as_pathway_db = TRUE)
```

The KidsFirst default can also load an external already-built GMT:

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

The current KidsFirst default is built from:

- `HOPE_pathway_database_without_KEGG_MEDICUS`: HOPE pathway database after
  removing `KEGG_MEDICUS*` pathways.
- `MSigDB_c2_cp_kegg_v7_canonical`: canonical KEGG pathways from MSigDB C2
  canonical pathways, added so standard names such as
  `KEGG_MAPK_SIGNALING_PATHWAY` are available.

With the current source table this corresponds to 8,785 HOPE-derived pathways
and 184 canonical MSigDB KEGG pathways after KidsFirst filtering.

For methods sections or audit trails, `kfirst_gosets_provenance()` returns a
citation-ready provenance statement, the component table, and the size filters:

```r
prov <- kfirst_gosets_provenance()

cat(prov$text)
prov$components
prov$filters
```

or rebuild the current KidsFirst convention from source GMTs:

```r
kfirst_db <- build_kfirst_default_pathway_db(
  hope_gmt = "path/to/pathway_database_HOPE.gmt",
  canonical_kegg_gmt = "path/to/c2.cp.kegg.v7.0.symbols.gmt",
  universe = interrogated_genes
)
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
