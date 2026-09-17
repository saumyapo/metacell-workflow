# Documentation

The code in `scripts/` is the pipeline. These pages are the reasoning behind
it: what a metacell is, why each parameter is set where it is, and how to tell
a finished run from a correct one.

Nothing here assumes prior experience with metacells. It does assume basic
familiarity with scRNA-seq (counts, UMIs, cell types, clustering). If that is
also new, [Background](01_background.md) lists the places to start.

## Reading order

Starting from nothing, read 01 to 03 in order, then run the pipeline and read
the rest as each step comes up.

| Page | Covers |
|:--|:--|
| [01 Background](01_background.md) | Why coarse-graining exists, where metacells sit among the alternatives, vocabulary |
| [02 What a metacell is](02_metacells.md) | The idea, how MC2 builds one, what it costs |
| [03 Preprocessing](03_preprocessing.md) | Gene deduplication, doublets, why V(D)J catches ones expression cannot |
| [04 Building metacells](04_building_metacells.md) | Cell QC, feature selection, divide and conquer, outliers |
| [05 Gene lists](05_gene_lists.md) | Excluded, lateral and noisy: three lists, three jobs |
| [06 Choosing parameters](06_parameters.md) | Graining level, UMI targets, pile size, module count |
| [07 Annotation](07_annotation.md) | Gene modules, blocks, marker scoring, manual curation |
| [08 Two passes](08_two_passes.md) | Why subset, the purity gate, stale field stripping |
| [09 Composition](09_composition.md) | Cell frequencies, why proportions need CLR |
| [10 Differential expression](10_differential_expression.md) | Pseudobulk, the four contrast families, the interaction problem |
| [11 Pathway enrichment](11_gsea.md) | Ranking statistics and why not fold change |
| [12 Clonal repertoire](12_clonal.md) | Clone calling, size scales, diversity, rarefaction |
| [13 Reading the diagnostics](13_diagnostics.md) | Every printed check, what it means, what to do |
| [14 Reproducibility](14_reproducibility.md) | Environments, containers, seeds, versions |
| [15 Outputs](15_outputs.md) | What each step writes and which step reads it |
| [16 Resources](16_resources.md) | External papers, tutorials and package docs |

## Shortcuts

- Never used metacells: [02](02_metacells.md). Nothing downstream makes sense first.
- Used MC2 before, want to know what this pipeline does differently:
  [05](05_gene_lists.md) and [06](06_parameters.md), where most of the decisions live.
- Run finished but the output looks wrong: [13](13_diagnostics.md).
- Looking for a file a script expects and cannot find: [15](15_outputs.md).
