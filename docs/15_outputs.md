# Outputs

What each step writes, and which step reads it. `<slug>` is `global` or the
subset slug; `<prefix>` is `FILE_PREFIX`, set in each script's config chunk.

Everything lands under `results/`, created on demand:

```
results/
  h5ad/<slug>/metacell/     AnnData objects, written by 01, 02, 04
  rds/<slug>/metacell/      Seurat objects, written by 03, 05, 10
  csvs/<slug>/metacell/     tables and workbooks
  pdfs/<slug>/metacell/     figures
```

## Objects

| File | Written by | Read by |
|:--|:--|:--|
| `project.merged_samples.h5ad` | 01 | 02 |
| `project.clean_cells.h5ad` | 02 | 03, 04 |
| `project.metacells.h5ad` | 02 | 03 |
| `project.subset.clean_cells.h5ad` | 04 | 05 |
| `project.subset.metacells.h5ad` | 04 | 05 |
| `<prefix>_MC_annotated.RDS` | 03, 05 | 06, 07, 08, 10, 12 |
| `<prefix>_MC_TCR.RDS` | 10 | 11 |

`clean_cells.h5ad` carries the per-cell metacell assignment and the
`lateral_gene` / `noisy_gene` / `rare_gene` masks. Both it and the metacell
object are needed by the annotation step.

## Tables

| File | Written by | Contents |
|:--|:--|:--|
| `project_<slug>_metacell_QC.xlsx` | 01, appended by 02 and 04 | Per-sample QC, UMI percentiles, attrition, outlier rate |
| `<prefix>_metacell_block_annotations.csv` | 03, 05 | One row per metacell: block, label, and the state columns |
| `<prefix>_module_block_map.xlsx` | 03, 05 | Module membership, block means, block to cell type mapping and scores |
| `<prefix>_per_sample_composition.xlsx` | 06 | Cell counts and percentages per sample and label |
| `<prefix>_composition_summaries.xlsx` | 06 | Donor-weighted means behind each bar figure |
| `<prefix>_CLR_timepoint.xlsx` | 07 | CLR table, model fits, contrasts |
| `<prefix>_pseudobulk_DEGs_<family>.xlsx` | 08 | One sheet per cluster |
| `<prefix>_DEG_run_manifest.xlsx` | 08 | Every attempted contrast, units per arm, status |
| `<prefix>_pseudobulk_GSEA_<family>.xlsx` | 09 | One workbook per family, unfiltered |
| `<prefix>_pseudobulk_GSEA_<family>_filtered.xlsx` | 09 | `\|NES\| > 1.5`, `q < 0.05` |
| `<prefix>_GSEA_run_manifest.xlsx` | 09 | One row per DEG csv |
| `<prefix>_clonal_TCR.xlsx` | 10 | Every descriptive repertoire table |
| `<prefix>_clonal_TCR_stats.xlsx` | 11 | Diversity, tracking, sharing, test families |

The two `block_annotations.csv` files are the handoff between the annotation
step and everything else. Step 04 reads the global one to subset on
`broad_label`; the composition, DE and clonal steps read labels off the
annotated RDS, which carries the same calls.

## Per-contrast DEG csvs

Step 08 also writes one csv per contrast, in a family folder, using the
filename grammar step 09 depends on:

```
<prefix>__<unit>__<contrast>.csv          ->  <prefix>__<unit>__<contrast>__GSEA.xlsx
```

One DEG csv gives one GSEA workbook with the same stem. A mismatch in counts
between the two manifests is the first thing to check when a result is missing.

## Figures

Named `<prefix>_<what>.pdf` in `results/pdfs/<slug>/metacell/`: the module
heatmap and annotation heatmap from 03 and 05, stacked and dodged composition
bars from 06, forest plots from 07, volcanoes from 08, top-pathway pages from
09, and the repertoire figures from 10 and 11.

Knitted HTML from each Rmd is worth keeping next to these. It records the
package versions and every diagnostic the run printed.
