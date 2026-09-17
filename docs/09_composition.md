# Composition

Steps 06 and 07. How much of each cell type is there, and does it differ.

## Count cells, not metacells

Every frequency in this pipeline is over cells, each carrying its metacell's
block label.

Counting metacells would be wrong twice. MC2 sizes metacells to a target cell
count, so their number is not proportional to abundance. And a metacell can
pool cells from several samples, so it has no single sample of origin.

The metacell is the annotation device. It is not the counting unit.

## Dropped labels

`Outlier`, `Unassigned` and `Undefined` are not cell types and are removed
before any proportion is computed.

That is only harmless if the drop rate is flat across samples. If one condition
loses many more cells than another, every proportion is shifted by that alone.
Step 06 prints the drop rate per sample and per condition first, and it is the
table to read before trusting any plot below it.

## `complete()` is load-bearing

```r
dplyr::count(sample_origin, cell_label) %>%
  tidyr::complete(sample_origin, cell_label, fill = list(n = 0))
```

A label absent from a sample is a real zero for that sample, not a missing row.
Without `complete()`:

- the mean for a rare label is taken over only the donors that happen to carry
  it, which inflates it
- per-group percentages stop summing to 100
- `position = "fill"` silently renormalises the bars so it cannot be seen

The scripts assert on this rather than trusting it:

```r
stopifnot(all(abs(tapply(d$mean_pct, d$grp, sum) - 100) < 1e-6))
```

Bars are drawn with `position = "stack"`, not `"fill"`, for the same reason. If
the means ever stop closing to 100, the bars should look wrong.

## Donor weighting

Bars are the mean across donors of each donor's own percentage, so every donor
carries equal weight regardless of how many cells were captured. Pooling cells
across donors instead lets whichever sample sequenced deepest set the answer.

## Why CLR

Proportions are compositional: constrained to sum to 100 within a sample. A
real increase in one population forces an apparent decrease in every other,
whether or not anything happened to those others. Testing raw proportions with
a t-test generates significant results for populations that did not change.

The centred log-ratio maps each sample out of the simplex into unconstrained
Euclidean space:

```
clr(x_i) = log( x_i / geometric_mean(x) )
```

After that, ordinary linear models are valid.

The interpretation changes, and this belongs in the figure legend. A CLR value
is abundance relative to the average population in the same sample, not an
absolute percentage. A positive CLR difference means that population grew
relative to the sample's typical population. It does not mean its absolute
count rose.

Alternatives worth knowing about: `propeller`
([Phipson et al. 2022](https://doi.org/10.1093/bioinformatics/btac582)) and
scCODA ([Büttner et al. 2021](https://doi.org/10.1038/s41467-021-27150-6)) both
handle compositional cell-type data with more machinery. CLR plus a mixed model
is the lighter option and is what this pipeline uses.

## Zeros

`log(0)` is `-Inf` and would poison every model that sample enters.

Zeros are replaced multiplicatively with half the smallest observed non-zero
proportion in the same sample, which is the standard compositional treatment,
and the count of replacements is printed so the scale of the correction stays
visible. The row is then reclosed to 100. CLR is invariant to closure so this
changes nothing numerically; it keeps the input table readable.

If the printed zero count is a large fraction of the table, the label set is
too fine for the data, and the honest fix is to merge labels rather than lean
harder on the imputation.

## The models

Which model depends on what varies where:

| Contrast | Model | Why |
|:--|:--|:--|
| Timepoints within a condition | `lmer(clr ~ timepoint + (1 \| patient_id))` | The same patient is sampled repeatedly, so patient is a random intercept |
| Conditions at one timepoint | `lm(clr ~ condition)` | One sample per patient at that timepoint, nothing to block on |

The mixed model falls back to `lm` when there is no repeated sampling to
exploit, and which was used is recorded in a `model` column.

> [!WARNING]
> `lmerTest` must be attached, not just `lme4`. Without it, `anova()` on an
> `lmer` fit has no p-value column, and `emmeans` falls back to asymptotic z
> inference, which is anti-conservative at n = 3 to 5.

Timepoints with fewer than `MIN_PATIENTS_PER_TP` patients are dropped before
fitting rather than filtered after. Dropped after, they still inflate the
residual degrees of freedom and produce contrasts leaning on a single sample.

Multiple testing is BH-corrected over the contrasts within each cluster.

## One `as.data.frame()` trap

```r
clr_long <- tibble::as_tibble(clr_tab, .name_repair = "minimal")
```

`as.data.frame()` on a matrix runs `make.names()` on the columns, which turns a
label like `CD4/CD8-amb Exhausted/Act` into `CD4.CD8.amb.Exhausted.Act` and
silently breaks every join back to the label set. `as_tibble` with minimal
repair keeps the names verbatim.
