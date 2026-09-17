# Building metacells

Steps 02 and 04. This page covers the mechanics common to both. What is
specific to the second pass is in [08](08_two_passes.md).

## Order of operations

Per-cell QC metrics are computed before any gene is dropped. This is not
cosmetic: percent mito needs the mitochondrial genes still present, and
`exclude_genes` removes them. Computed afterwards,
`str.startswith("MT-")` matches nothing and returns a clean 0 percent for every
cell.

```
compute total_umis, n_genes, percent_mito
  -> exclude_genes
  -> exclude_cells
  -> extract_clean_data
  -> mark lateral / noisy
  -> select features
  -> divide and conquer
  -> collect
```

## Mitochondrial fraction

```python
mito_mask = merged.var_names.str.startswith("MT-")
```

Not `MT-` or `MTRNR`. The `MTRNR2L*` genes are nuclear-encoded humanin-like
pseudogenes, not mitochondrial, and including them inflates the percentage.
They belong on the excluded list for a different reason (ambient,
stress-associated), which is where this pipeline puts them.

## Cell filtering

Two things act.

**A UMI floor.** `min_umi_floor: 500`, a safety net rather than the main
filter. The upper bound is set deliberately high, because a high-count cell is
better handled as a doublet question than as a QC threshold.

**Excluded gene fraction.** A cell whose mitochondrial, haemoglobin and
`MALAT1` counts exceed `max_excluded_gene_fraction` of its total is dropped.
This filter does the real work. It is composite by design: it catches dying
cells, ambient-dominated droplets and red cell contamination with one
threshold, and it scales with library size instead of being a fixed cutoff.

### MAD and miQC

Both scripts contain implementations of MAD-based per-sample outlier calling
and a Python reimplementation of
[miQC](https://doi.org/10.1371/journal.pcbi.1009290), currently commented out.

They are kept rather than deleted because they are a real alternative, not dead
code. MAD thresholds
([Heumos et al. 2023](https://doi.org/10.1038/s41576-023-00586-w)) are
data-driven rather than copied cutoffs, and miQC accounts for cell complexity,
so a high mito fraction in a low-complexity cell is flagged while the same
fraction in a high-complexity cell is not.

They are off because the composite excluded-fraction filter already removes the
same cells in this data, and stacking three filters makes the attrition table
hard to reason about. Switching them on OR-s the masks into `excluded_cell`,
and the attrition table has columns waiting for them.

> [!TIP]
> Use MAD with `nmads=5` on raw MAD units, matching the sc-best-practices
> convention, or `nmads=3` on SD-equivalent units, matching scater. Do not
> combine `normal_scale=True` with `nmads=5` unless 5 SD-equivalents is meant,
> which is a band roughly 7.4 raw MADs wide that no standard uses.

## Feature selection

```python
selected = mc.pl.extract_selected_data(
    adata=clean, min_gene_relative_variance=None, random_seed=RANDOM_SEED)
```

`min_gene_relative_variance=None` lets MC2 pick the cutoff from the data rather
than imposing one. The leak check that follows is the important part and is
covered in [05](05_gene_lists.md).

## Divide and conquer

```python
mc.pl.compute_target_pile_size(...)
mc.pl.divide_and_conquer_pipeline(...)
```

`compute_target_pile_size()` has to run first: it writes the pile size onto the
object and the pipeline reads it. The parameters go to both because both stages
use them.

Wrap the pipeline in `mc.ut.progress_bar()`. On a few hundred thousand cells it
runs for hours, and the bar is the only signal that it is still working.

Set the worker count from the scheduler rather than hardcoding it:

```python
n_processors = int(os.environ.get("LSB_DJOB_NUMPROC", os.cpu_count() or 4))
mc.ut.set_processors_count(n_processors)
```

Both container definitions also pin `OMP_NUM_THREADS=1`, so BLAS does not
oversubscribe on top of MC2's own parallelism.

## Outliers

Cells MC2 could not place get `metacell = -1`.

```
Metacells: 6367 | outlier cells: 12841 (5.05%)
```

5 to 10 percent is normal. The rate usually rises on a subset, because the
population is more homogeneous and there is less signal to separate on. Above
roughly 10 percent on a subset suggests `target_metacell_size` is too large for
how much structure is left.

Outliers are reported per sample in the attrition table. One sample with a much
higher outlier rate than the others is a batch signal worth chasing, not a
per-cell problem.

## Rare gene modules

MC2 detects modules of genes expressed in a very small fraction of cells at
high level, which is how genuinely rare populations survive a process that
would otherwise absorb them.

Both scripts pull the lists out and write them to the QC workbook, for two
reasons. Rare modules identify rare cell types, and they flag genes that should
probably have been lateral: a rare module that turns out to be immunoglobulin
or a stress programme means the lateral list has a gap.

> [!WARNING]
> The rare gene detector has a cap of 500 candidate genes. If contaminant genes
> fill it, the genes that would have separated a real rare population never get
> considered. Step 04 checks whether the cap is being hit:
>
> ```python
> candidate = (n_expressing / n < 0.001) & (gene_max >= 7)
> print(f"candidate rare genes: {int(candidate.sum())} (cap is 500)")
> ```
>
> Under the cap, this is not a concern.

Step 04 also cross-tabulates rare modules against global block and label, so a
module driven by contamination from one block can be traced and, if needed,
those cells dropped before the subset is rebuilt.

On a subset an empty rare-module result is common and informative: the rare
populations the global run was detecting are no longer in the object.

## What gets written

| File | Contents |
|:--|:--|
| `*.clean_cells.h5ad` | Per-cell object with metacell assignment and the `lateral_gene` / `noisy_gene` / `rare_gene` masks. Steps 03 and 05 read the masks. |
| `*.metacells.h5ad` | The collected metacell object, one row per metacell. |

Both are needed downstream. The masks on the cell object are what let step 03
exclude lateral genes from the module HVG selection.
