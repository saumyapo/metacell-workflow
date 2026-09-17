# Reading the diagnostics

Every check the scripts print, what it means, and what to do when it looks
wrong. This is the page for a run that finished and looks off.

## Step 01

| Output | Healthy | If not |
|:--|:--|:--|
| Duplicated genes after merge | Empty vector per sample | The merge did not take. Stop; every total downstream is wrong. |
| TCR recovery percent | 75 to 90 | Below that on one sample, check the V(D)J library resolved to the right folder. Below that everywhere, check the chemistry. |
| Dual-beta raw vs supported | Supported well below raw | Equal means the ambient filter is not firing; check `min_second_beta_umis`. |
| TCR doublets removed | Well under the 10x expectation | Above it means the filter is over-calling. Raise `min_second_beta_frac`. |
| Retention percent per sample | Flat across samples | One sample far below the others is a batch problem, and every proportion downstream inherits it. |

## Steps 02 and 04

**Cells and genes retained.**

```
Cells: 254331 -> 251570 (98.91%)
Genes: 33344 -> 33331 (99.96%)
```

Cell retention in the high 90s is normal. A large drop means the
excluded-fraction threshold is too tight, or a sample is genuinely poor.

**Attrition table.** Per sample, per filter stage. Read across rows: one sample
losing far more than the others at one stage names the problem.

**Selected feature genes.** A few thousand globally, fewer on a subset. Below
about 200 on a subset the module step will not be meaningful; the cause is
usually an over-aggressive purity gate or a lateral list that ate the
informative genes.

**Lateral leak check.**

```
Lateral genes that leaked into selection (should be empty): []
```

Anything printed here means the lateral marking did not take effect, usually
applied to the wrong object or after selection. Fix and rerun; do not proceed.

**TCR V/J check** (subset only). Anything printed here on a lineage subset
means clonal signal is driving the similarity graph. Stop.

**Outlier rate.**

```
Metacells: 6367 | outlier cells: 12841 (5.05%)
```

5 to 10 percent normal, rising on a subset expected. Above 10 percent on a
subset means `target_metacell_size` is too large for the remaining structure.

**Rare gene cap.**

```
candidate rare genes: 312 (cap is 500)
```

Under the cap, no concern. At the cap, contaminant genes may be crowding out
the genes that would separate a real rare population.

## Step 04, the purity gate

**Contaminant distribution within labelled blocks.** Printed at several
percentiles. Set `purity_other_max` off the knee rather than at a round number.

**Per-block failure rate.** Concentrated failure, one block at 80 percent,
means that block's label is wrong and belongs fixed in step 03. Diffuse failure
across blocks is stray metacells, which is what the gate is for.

**Mean panel of dropped metacells.** They should look foreign, high on a
contaminant panel. Low on everything including the lineage marker means the
floor is too high and shallow real cells are being lost.

**Cells in flagged blocks.**

```
cells in flagged blocks: 1,204 (0.47% of clean cells)
```

At a fraction of a percent, the ambiguity question is academic and the default
subset is fine as it stands.

## Steps 03 and 05

**`range(x)`.** Read before setting the HVG mean threshold. It shifts between
the global object and every subset.

**Threshold sweep.** Pick from the curve. Target 1,500 to 2,500 module genes.

**Module size distribution.**

```r
table(cut(table(res2_cut), breaks = c(0,1,3,5,10,20,Inf)))
```

Many 1-to-3-gene modules means `NCLUST` is too high. Drop to 60 to 80.

**Blocks with fewer than 10 metacells.** If several appear, lower `column_km`
and redraw the heatmap.

**`max(abs(gz)) < 10`.** Asserted. A z above 10 means a scale estimate
collapsed, which happens when a gene is on in only one or two blocks.

**Metacell accounting.**

```
Metacells in module heatmap (thmz):  6367
Metacells assigned to >=1 cell:      6367
Cells total:                       251570
  in metacell:                     238729
  MC2 outliers:                     12841
  Unassigned (MC not in heatmap):       0
```

`Unassigned` should be 0. Anything else means metacells were dropped between
the metacell object and the heatmap matrix, and the annotation is incomplete in
a way nothing downstream will flag.

## Steps 06 and 07

**Drop QC table.** Read before any plot. Drop rate should be flat across
conditions.

**Percentages closing to 100.** Asserted rather than checked. If it fires,
`complete()` did not do its job.

**Structural zeros.** A large count relative to the table means the label set
is too fine, and the honest fix is merging labels.

**CLR centring.** Asserted: every sample's CLR values sum to zero.

**`model` column.** Records `lmer` or `lm` per cluster. A cluster falling back
to `lm` where the mixed model was expected means the repeated sampling was not
there for that cluster.

## Step 08

**Run manifest.** The most useful sheet in the workbook. Check:

- `skipped` rows and the reason
- `exploratory` rows, which rest on fewer than 3 units per arm
- `n_test` and `n_ref` matching the intended design

**`n_genes_tested`.** Should be in the thousands. Low hundreds means
`filterByExpr` removed almost everything, usually because the pseudobulk units
were built from too few cells.

## Step 09

**GSEA manifest.** One row per DEG csv. Statuses:

| Status | Meaning |
|:--|:--|
| `ok` | Tested |
| `skipped` | Ranked list shorter than `MIN_GENES` |
| `not rankable` | Omnibus F table, expected |
| `no result` | Ran, returned nothing |
| `no gene column` | Malformed input |

Row count should equal the number of DEG csvs. A mismatch is the first thing to
check when a result is missing.

## Steps 10 and 11

**V(D)J recovery per sample.** 75 to 90 percent.

**`clonalFrequency` is per-sample.** Prints `TRUE`. If `FALSE`, the `group.by`
argument did not take and every expansion statement is at the wrong level.

**Contig barcodes before and after filtering.** A sample dropping to zero means
the barcode reformat is wrong, not that V(D)J failed.

**Samples below `MIN_VDJ_CELLS`.** Reported, not silently dropped. If most
samples fall below it, lower it and accept wider intervals rather than
analysing three samples.

**Chain dropout percent.** Flagged above `DROPOUT_FLAG_PCT`. A high
alpha-dropout sample will bias any paired-chain clone call.
