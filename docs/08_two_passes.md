# Two passes

Metacells are built twice: once on everything, once inside one lineage. This
page is why, and what has to be handled at the boundary.

## Why not once

A single global run at fine grain sounds simpler and works badly. Feature
selection is global. With monocytes and T cells in the same object, the genes
that separate them dominate the variance, so those are the genes that get
selected and those are the genes the similarity graph runs on. Within-T-cell
structure (naive, memory, effector) is a smaller signal sitting underneath and
does not make the cut.

Subsetting removes the between-lineage variance. Feature selection re-runs on
what is left and picks up the state-level genes. Far fewer selected genes than
globally is the point, not a failure.

The same structure repeats for as many lineages as there are. Each subset is
its own run of 04 and 05 with its own config block.

## The purity gate

Block labels apply to every metacell in a block. If a block of 600 metacells is
called "T cells", a stray monocyte metacell inside it inherits that label and
walks into the subset unchallenged. That is what puts monocyte, B and cDC gene
modules into a T cell subset's module run, and by then it is hard to trace
back.

So step 04 scores every metacell on a small panel per lineage and gates on
absolute values:

```python
T_MIN, OTHER_MAX = 0.3, 0.4
mc_pure = (panel_mc["T"] >= T_MIN) & (other <= OTHER_MAX)
```

**Absolute ceilings, not a margin.** A T-minus-other margin punishes shallow
metacells: T = 0.40 with Mono = 0.07 is clean but has a gap of only 0.33, and
that gap is a depth artefact rather than an identity problem. Two independent
thresholds separate "has lineage signal" from "has foreign signal", which are
different questions.

Three diagnostics decide whether the gate is set right:

**Per-block failure rate.** Concentrated failures, one block failing at 80
percent, mean a block label is wrong and should be fixed upstream rather than
gated away. Diffuse failures across many blocks are stray metacells, which is
what the gate exists for.

**Mean panel of dropped metacells.** Dropped metacells should look foreign,
high on Mono or B. Low on everything including the lineage marker means the
floor is too high and real, shallow cells are being lost.

**The contaminant distribution within labelled blocks.** Printed at several
percentiles, so `OTHER_MAX` can be set off the knee rather than picked.

> [!TIP]
> Do not relabel a borderline block to move more cells into the subset on the
> theory that more cells is better. A contaminated block hurts the DE more than
> a slightly undersized subset does.

### One panel deliberately not gated on

NK markers are reported but never gated. `KLRD1` is CD94 and sits on CD8
effector-memory, Temra and gamma-delta T cells, so gating on it rejects real T
blocks. The T versus NK call is made at block level in step 03 by an explicit
CD3 veto, which is the right place for it.

## Stripping stale MC2 fields

MC2 does not clear its own annotations. If `metacell`, `metacell_name`,
`excluded_gene`, `lateral_gene`, `selected_gene`, `rare_gene_module` and the
rest survive from the global run into the subset object, `exclude_genes` and
the divide-and-conquer pipeline will reuse the global masks instead of
recomputing them.

The failure mode is bad: the run completes, the output looks structurally
normal, and the metacells were built on global feature selection. Nothing flags
it.

Step 04 renames anything worth keeping with a `global_` prefix, then deletes
the originals along with `obsp`, `varp` and `uns`:

```python
keep_as_global = ["metacell", "metacell_name", "metacell_level",
                  "dissolved", "most_similar", "most_similar_name"]
```

`global_broad_label`, `global_mc_block` and `global_anno_margin` are kept
deliberately, so any cell in the subset can be traced back to its parent
metacell in the global run.

## Metacell name collisions

Global metacells are `M1`, `M2`, ... Subset metacells are also `M1`, `M2`, ...
Everything downstream keys on `metacell_name`, so joining or `rbind`-ing a
subset table with a global one would silently match the wrong rows.

Step 04 prefixes at the source:

```python
metacells.obs_names = [f"{MC_PREFIX}{n}" for n in metacells.obs_names]
```

Prefixing once, where the names are created, means every lookup table, Seurat
object and CSV downstream carries unambiguous IDs without anything else having
to know about it. Set `mc_prefix` per lineage in `config/config.yaml`.

## No QC on the subset

Step 04 runs no per-cell QC. Cells arriving there already passed the global
filters, and re-filtering on a subset would recompute MAD thresholds against a
different population and silently drop real biology. A naive T cell has a
genuinely lower gene count than a monocyte; against the global distribution
that is fine, against a T-only distribution it starts looking like an outlier.

`excluded_cell` is set to all-`False` explicitly, so `extract_clean_data()` has
the field it expects and the gene filter is the only thing acting.

The one filter that is added is zero-count genes. Subsetting to one lineage
leaves several thousand genes at exactly zero UMIs, and those destabilise the
variance-mean LOESS fit in step 05.

## What changes in step 05

Same machinery as step 03, different tuning:

| | Global (03) | Subset (05) |
|:--|:--|:--|
| HVG mean threshold | -9 | Re-derived from `range(x)` |
| HVG varmean threshold | 3.67 | ~1.3, lower ceiling on one lineage |
| Modules (`NCLUST`) | 100 | 90 |
| Blocks (`column_km`) | 40 | 25 |
| Marker panel | Broad lineage | Lineage states |
| Scoring centre | Mean across blocks | Median across blocks |

That last row matters. Blocks are not balanced across states: with 10 naive
blocks out of 25, a mean across blocks sits inside the naive population and
naive blocks score about zero on their own markers. The median is not dragged
by the dominant state. The scale stays SD, because MAD collapses toward zero
for genes on in only one or two blocks and sends their z into the hundreds.
