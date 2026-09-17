# Annotation

Steps 03 and 05. Metacells go in, cell type labels come out. The route is gene
modules, then metacell blocks, then marker scoring on blocks.

## Why modules rather than markers directly

Scoring each metacell on a marker panel and taking the argmax works badly. A
single metacell is still noisy enough that one dropped marker flips the call,
and a marker panel only finds what someone already thought to look for.

Modules are found from the data. Genes that co-vary across metacells are
grouped, each metacell is scored on every module, and the blocks that emerge
are defined by which programmes they run. The marker panel is then applied to
blocks, which are averages over tens to hundreds of metacells and far more
stable. It also means an unexpected programme shows up as a module with no
label, which is a lead rather than a silent miss.

## Selecting module genes

Genes are placed on a mean-dispersion plot: log10 mean expression against
log2(variance / mean).

Dispersion rather than variance, because variance rises with the mean, so
ranking on variance alone just selects the most abundant genes.

A LOESS curve is fitted to the minimum dispersion in each expression bin. That
minimum represents technical noise at that abundance: the dispersion of a gene
that is not variable. Genes are selected when they sit a set distance above the
curve and above a floor on mean expression.

> [!WARNING]
> MC2 normalises each metacell to sum to 1, so mean expression is a mean of
> frequencies and log10 of it is negative. The range of that value shifts
> between the global object and a subset, because frequencies renormalise
> within whatever object is in hand. Thresholds do not transfer. Both scripts
> print `range(x)` and include a threshold sweep for this reason.

Lateral and noisy genes are removed from the selection afterwards. Expect a
larger proportional drop on a subset: TCR V/J and ribosomal genes make up more
of the top-dispersion set once between-lineage signal is gone.

## Modules

Pearson correlation between all selected gene pairs, Euclidean distance on the
correlation matrix, complete-linkage hierarchical clustering, `cutree` at
`NCLUST`.

Distance on the correlation matrix rather than on expression directly means two
genes are close when they correlate similarly with every other gene, not just
with each other. That picks up genes in the same programme even when they are
not individually well correlated.

Each metacell's module score is the fraction of its total log-expression
attributable to that module's genes, then z-scored across metacells.

## Blocks

The module-score matrix, modules by metacells, is drawn as a ComplexHeatmap
with `column_km` k-means groups. Those groups are the blocks.

`column_km` is the resolution knob, tuned by looking at the heatmap and at the
metacells-per-block table. The scripts warn for blocks holding fewer than 10
metacells, since those give unstable lineage scores.

`blockmeans` gives module by block mean z, and the top five modules per block
are printed. This is the table to read when a block's call looks wrong: it says
which programmes the block is actually running.

## Scoring blocks

Two gene sets per lineage:

- **`panel`** is what gets displayed on the annotation heatmap. Generous,
  includes markers that appear under more than one cell type.
- **`sig`** is what gets scored. Every gene appears in exactly one signature.

The separation matters. Anything a reference like Azimuth PBMC L2 lists under
two different subsets (`IL7R`, `KLRB1`, `NKG7`, `CCL5`, `KLRD1`, `CTLA4`,
`EOMES`) is display-only. A gene that cannot separate two states in the
reference will not separate them here, and putting it in `sig` adds noise to
both.

Genes are z-scored per gene across blocks, then averaged within a signature.
Averaging raw log-CP10K instead lets one high-expression gene carry a whole
signature, which is how `HLA-DRA` alone ended up defining an activation
signature in an early version.

The call is the argmax, subject to:

- `best < BEST_MIN` becomes `Undefined`. Nothing fits.
- `margin < MARGIN_MIN` becomes `Undefined`. Two unrelated states fit equally.
- Adjacent pairs are exempt from the margin rule. States that share a programme
  legitimately tie, and a near-tie between them is overlap rather than
  ambiguity. Without this, gamma-delta and NKT blocks land in `Undefined`
  against a cytotoxic signature they genuinely resemble.

### The CD3 veto

The T versus NK boundary is decided on absolute expression, not on z scores, in
both directions:

```r
CD3_MIN <- 1.0
bad_T  <- call == "T cells" & (cd3 < CD3_MIN | nk > cd3)
bad_NK <- call == "NK"      & cd3 >= CD3_MIN & cd3 > nk
```

The two lineages share too much of a cytotoxic programme for a relative score
to settle it, and CD3 is the definitional marker. Doing this at block level
means the subset step downstream does not have to re-litigate it.

### Orthogonal axes

CD4 versus CD8 is scored separately from cell state and scaled per column,
because the two are only compared against each other. Folding it into the main
`max.col()` would force every CD8 block to be labelled CD8 regardless of
whether it is naive, GZMK+ or Temra.

## Manual curation

Automatic calls are sometimes wrong, and the pipeline is built so they can be
fixed without losing the audit trail.

```r
BLOCK_RELABEL <- c()   # e.g. c("7" = "T cells", "5" = "Granulocytes")
```

`block_anno` always holds the automatic call. Curation writes
`block_anno_final`, and a `label_source` column records which was used. Nothing
overwrites the automatic result, so a reviewer can see both.

> [!WARNING]
> Do not relabel `Undefined` blocks wholesale. If nothing scored above zero,
> the honest read is that the block has no confident identity.

Two things justify a relabel: the top modules for that block clearly name a
programme the signature missed, or a definitional marker is present and the
score failed for a technical reason such as dropout in resting cells. Both are
recorded as comments next to the entry, because a hand-curated label is a
methods decision and has to be reportable.
