# Choosing parameters

Read this before changing anything in `config/config.yaml`.

## Graining level: `target_metacell_size`

The graining level is gamma = cells / metacells, set directly by
`target_metacell_size`. Everything else about the run follows from it.

[Bilous et al. 2022](https://doi.org/10.1186/s12859-022-04861-1) benchmarked
this and found gamma between 10 and 50 keeps results faithful: in that range,
over 75 percent of single-cell-level differentially expressed genes were still
recovered from metacells. Above 50 they stopped testing, which is a reason to
treat 50 as a ceiling rather than evidence that 100 is fine.

Within that range:

| Smaller gamma (20 to 30) | Larger gamma (40 to 50) |
|:--|:--|
| More metacells, finer structure | Fewer metacells, deeper each |
| Rare states survive as their own units | Rare states get absorbed |
| Noisier per-metacell profiles | Cleaner profiles, better module scores |
| Slower | Faster |

Gamma depends on cell counts, not cell type. What changes between lineages is
depth, not how many cells need pooling.

## UMI target: `target_metacell_umis`

MC2 takes a cell target and a UMI target, and whichever binds first is the one
that acts. So the two have to be set together or one of them does nothing.

```
target_metacell_umis  ~=  target_metacell_size x median UMIs per cell
```

Step 04 prints exactly this, on the post-exclusion totals for the subset:

```
n_cells=120,317  median=3,120  IQR=1,980-4,870
  size=24 -> umis=75,000
  size=36 -> umis=112,000
  size=48 -> umis=150,000
```

Set the UMI target from that table, not from the global run's value.

**UMI target far above size x median**: the size target binds and metacells
come out roughly the requested size. This is usually what is wanted.

**Far below**: the UMI cap binds and metacells come out much smaller than the
size target, which is confusing because the size was set and something else
came back.

> [!WARNING]
> `total_umis` on the subset object is the post-exclusion total. The global
> `exclude_genes` already removed mitochondrial, haemoglobin and `MALAT1`
> counts, so it is meaningfully lower than raw library size. Carrying the
> global UMI target onto a subset without recomputing is a common way to end up
> with the wrong binding constraint.

## Depth differs by lineage

Same gamma, different UMI target, because library depth is not uniform:

| Lineage | Consideration | Suggested size | UMIs |
|:--|:--|:--|:--|
| T cells | Usually the largest subset, gamma can sit high | 40 to 50 | size x median |
| Myeloid | Deeper libraries, fewer cells | 30 to 40 | size x median |
| B cells | Usually the smallest of the three | 24 to 32 | size x median |

## Piles: `min_pile`, `max_pile`, `target_metacells_in_pile`

Piles are how divide and conquer chunks the data. They are a computational
device and the reconciliation step is meant to remove their trace, but the
sizes still matter at the edges.

`compute_target_pile_size()` derives the actual pile size from
`target_metacells_in_pile` and clamps it between the min and max. Leaving
`target_metacells_in_pile` at 100 is reasonable and rarely worth touching.

The one real decision is whether to use divide and conquer at all. Step 04
switches on cell count:

```python
use_dac = clean.n_obs >= 2 * MIN_PILE
```

Below two piles' worth of cells it calls `compute_metacells()` directly, which
is faster and avoids pile-boundary artefacts on a population that would
otherwise be split arbitrarily for no benefit. For small subsets (pDCs,
platelets, plasma cells) the direct call is the right one. Where fine
separation matters, lowering `min_pile` to force divide and conquer can give
more precise calls.

## Module count: `NCLUST`

The number of gene modules cut from the gene-gene correlation tree in steps 03
and 05. Default 100 globally, 90 on a subset.

The constraint that gets missed: gene-gene correlations here are computed
across metacells, so their stability scales with the number of metacells, not
the number of cells. Below roughly 200 metacells the correlations get noisy and
`NCLUST` should come down with them.

The diagnostic is the module size distribution, printed in step 05:

```r
table(cut(table(res2_cut), breaks = c(0,1,3,5,10,20,Inf)))
```

Many modules coming back with 1 to 3 genes means `NCLUST` is too high for this
subset. Drop to 60 to 80 and rerun.

## HVG thresholds: `inVarMean_MeanThresh`, `inVarMean_varmeanThresh`

These select the genes that go into module construction. The mean threshold is
a floor on log10 mean expression; the varmean threshold is how far above the
fitted LOESS noise curve a gene has to sit.

> [!WARNING]
> These must be re-tuned per subset. MC2 normalises each metacell's expression
> to sum to 1, so frequencies renormalise within whatever object is in hand.
> The range of log10(mean) shifts between the global object and a subset, and a
> threshold copied across will not select the same number of genes. Step 05
> prints `range(x)` for this reason. Read it first.

Both scripts include a sweep, so the threshold comes off the curve rather than
a guess:

```r
thresholds <- seq(0.5, 1.5, by = 0.05)
setNames(sapply(thresholds, n_genes), thresholds)
```

Aim for roughly 1,500 to 2,500 genes after lateral and noisy genes are removed.
Below about 200 the gene-gene correlation step stops being meaningful, and the
script asserts on it.

Expect substantially fewer selected genes on a subset than globally.
Between-lineage variance is gone and what remains is within-lineage structure,
which is a smaller signal. That is expected. If it drops below about 200 the
problem is usually upstream: too aggressive a purity gate, or a lateral list
that ate the informative genes.

## Number of blocks: `column_km`

How many groups the module-score heatmap is cut into. This is the closest thing
in the pipeline to a clustering resolution parameter, and like resolution it is
tuned by looking.

40 was chosen on roughly 6,000 global metacells. On a subset with a fraction of
that, 40 blocks leaves too few metacells per block for marker scoring to be
stable, so step 05 starts at 25.

The check to read:

```r
cat("blocks with <10 metacells:", ...)
```

If several appear, lower it and rerun the heatmap chunk.

## Summary

| Parameter | Default | Re-tune per subset | Read this first |
|:--|:--|:--|:--|
| `target_metacell_size` | 72 global, 24 subset | Yes | Cell count |
| `target_metacell_umis` | 240k / 160k | Always | Median post-exclusion UMIs |
| `target_metacells_in_pile` | 100 | Rarely | |
| `min_pile` / `max_pile` | 2000-16000 | Sometimes | Whether divide and conquer is wanted at all |
| `NCLUST` | 100 / 90 | Yes | Metacell count, module size table |
| `inVarMean_*` | 3.67 / 1.3 | Always | `range(x)`, threshold sweep |
| `column_km` | 40 / 25 | Yes | Metacells per block |
| `purity_min` / `purity_other_max` | 0.3 / 0.4 | Yes | Contaminant score distribution |
