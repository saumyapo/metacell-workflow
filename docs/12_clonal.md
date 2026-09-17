# Clonal repertoire

Steps 10 and 11. Step 10 attaches the repertoire and produces descriptive
figures. Step 11 does diversity, tracking and statistics from the object step
10 wrote.

The structure of step 10 follows the
[rnabio.org TCR workflow](https://rnabio.org/module-08-scrna/0008/06/01/Clonality/);
the analytical choices below are this pipeline's.

## Clone calling

```r
CLONE_CALL <- "strict"
```

V and J gene plus the nucleotide CDR3, on both chains. The most conservative
option scRepertoire offers: two cells share a clone only if their receptors are
identical at nucleotide level.

Looser calls (amino acid CDR3, or one chain only) merge convergent
recombination events into one clone. Within a patient that is usually the wrong
call. Across patients it is a different question, handled separately as public
clone detection in step 11, which is the only place `CTaa` is used.

## Filter contigs to cells in the object

Contig files contain barcodes that never reached the analysed object: cells
that failed RNA QC, called doublets, cells that landed in another lineage
subset.

Left in, every repertoire summary describes the raw contig file rather than the
cells being analysed, and the clone size denominators stop agreeing with the QC
table. Filtering before `combineExpression` rather than only inside it gives
one denominator everywhere.

`clonalFrequency` and `cloneSize` shift slightly compared to an unfiltered run.
That is the intended effect, not drift.

## V(D)J recovery is never 100 percent

Expect 75 to 90 percent. Cells with low TCR mRNA fail contig assembly, and
innate-like cells make no productive TCR reads at all. The per-sample recovery
table is printed so an outlier sample is visible; a sample far below the others
usually points at a V(D)J library problem rather than at biology.

`filterMulti = TRUE` keeps the highest-UMI contig per chain rather than
discarding the barcode, which is the right call after step 01's dual-TRB
doublet filter has already removed the real doublets. `removeNA = FALSE` keeps
single-chain barcodes, because they are still clones.

## Chain dropout

`CTstrict` concatenates both chains, so a cell that lost its alpha chain gets a
different string from its own sisters and is counted as a separate clone. That
inflates clone counts, deflates the largest clone sizes, biases clonality
downwards and under-detects persistence, at a rate that differs several-fold
between samples.

Step 10 measures it and carries `chain_status` and a beta-only clone id
(`clone_trb`) on the object, so step 11 can run the sensitivity analysis
without re-reading contigs. Orphans are not merged into their beta-matched
parent, because that would invent lineage calls the data does not support.

## Clone size: three scales

**Count bins** are absolute cell counts, tuned to PBMC where a 20-cell clone is
already large. Not comparable across samples of different depth: a 30-cell
clone is 3 percent of a 1,000-cell sample and 0.3 percent of a 10,000-cell
sample, and lands in the same bin either way.

**Occurrence bins** are a clone's share of its own sample's repertoire, as a
percentage. Depth-independent by construction, following
[Mhanna et al. 2024](https://doi.org/10.1016/j.crmeth.2024.100738) Fig 1C. This
is the scale to use for anything compared across samples.

**Expanded versus unique**, the `>= 2` cell dichotomy. The one clone-size
statement that involves no arbitrary cutoff at all.

Neither binned scale is derived from anything. They are display conventions,
which is why the rank-abundance figure is shown alongside them.

> [!WARNING]
> `group.by = "none"` computes clone frequency within each sample, which is
> what a timepoint-level expansion statement needs. The script verifies this
> rather than trusting the argument:
>
> ```r
> cat("clonalFrequency is per-sample:", all(.chk$n_in_sample == .chk$clonalFrequency))
> ```
>
> Patient-level frequency is computed separately, and it is the denominator for
> anything tracking a clone across visits.

## Barcode format

`combineTCR` writes `SAMPLE_BARCODE`. Step 01 wrote cell names as
`BARCODE_SAMPLE`. They have to be reconciled:

```r
d$barcode <- sub("^(.*?)_(.*)$", "\\2_\\1", d$barcode)
```

A sample losing every barcode after the join means this reformat is wrong, not
that V(D)J failed. The script checks for that specifically.

## Diversity needs rarefaction

Every diversity index is depth-dependent. Samples that differ several-fold in
TCR+ cell count cannot be compared on raw Shannon entropy, and no index is
immune to this.

Step 11 rarefies to a common cell count before computing anything:

```r
MIN_VDJ_CELLS <- 750   # below this, no diversity estimate
RAREFY_TO     <- 600   # NA means the smallest retained sample
```

Both are judgement calls. Raising `MIN_VDJ_CELLS` tightens the estimates and
drops samples; `RAREFY_TO` should sit at or below the smallest retained sample.
Anchoring it to the smallest sample automatically is convenient and lets one
shallow sample set the resolution for everything, so it is left explicit.

The rarefaction curves, observed and extrapolated, are the figure that
justifies this rather than a bare threshold.

Chao1 is not computed. On singleton-dominated TCR data it was unstable across
rarefactions, and sitting in the same BH family it inflated the adjusted
p-values of every other metric.

## The metric set

| Metric | Sensitive to |
|:--|:--|
| Renyi profile | The whole curve. Shannon is one point on it, at alpha = 1 |
| Shannon | The typical clone |
| Richness (alpha = 0) | How many distinct clones |
| Berger-Parker (alpha to Inf) | The single dominant clone |
| Gini | Inequality of the size distribution |

Plotting the whole Renyi profile says whether a difference is about richness,
about the typical clone, or about the dominant clone, and it stops the reader
arguing about which single index was chosen. Two repertoires whose curves cross
are not orderable by any single diversity index, which is worth knowing before
reporting one.

Gini is not redundant with normalised entropy. It responds differently to a
long tail and is the clonality measure used in several recent repertoire
papers.

## Tracking and sharing

**Across visits within a patient.** Morisita is the one to report: it is
abundance weighted, where Jaccard is dominated by singletons and can read 0.99
for two samples that share none of their dominant clones. A clone is persistent
when it is seen at two or more visits of the same patient. The clone-count
biplot shows the same information per clone: shared clones sit on the diagonal,
clones private to one visit lie along an axis.

**Between states.** Tversky asymmetric index, drawn as one square heatmap with
the diagonal masked. Asymmetry is informative: a small effector state can be
almost entirely contained in a large memory state without the reverse being
true, and a symmetric index averages that away.

## The replicate, again

Every metric is reduced to one number per sample, then tested with the patient
as the unit, using the same four families as step 08
([10](10_differential_expression.md)).

Some published analyses pool clones across patients for sharing and inequality
panels. This pipeline keeps the patient as the unit throughout: with a handful
of patients and uneven visit counts, pooling lets the deepest patient set the
answer.

With small arms the smallest attainable Wilcoxon p is bounded by the group
sizes alone, so the condition family is reported with effect sizes and an
explicit `powered` flag rather than as a test that can succeed.

## Rank-abundance

Rank against frequency on log-log axes. A steeper head means a more clonally
dominated repertoire, and it depends on no binning choice at all. When a binned
figure and the rank-abundance curve disagree, the curve is the one to trust.
