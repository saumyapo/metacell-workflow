# Differential expression

Step 08.

## Why pseudobulk

The replicate is the donor, not the cell. Cells from one donor are not
independent observations of a condition effect; they share that donor's
genotype, batch, dissociation and disease state.

Testing single cells directly treats several thousand correlated observations
as several thousand replicates. The effective sample size is inflated by orders
of magnitude and the p-values do not control the false positive rate at all,
not approximately, not conservatively. This is well established
([Squair et al. 2021](https://doi.org/10.1038/s41467-021-25960-2);
[Zimmerman et al. 2021](https://doi.org/10.1038/s41467-021-21038-1);
[Murphy and Skene 2022](https://doi.org/10.1038/s41467-022-35519-4)).

Counts are summed to one profile per replicate unit before any test runs.

## Which unit is the replicate

It depends on the contrast, and getting this wrong is the second most common
mistake after not pseudobulking at all.

| Contrast | Varies | Unit | Design |
|:--|:--|:--|:--|
| Case vs control, pooled | Between patients | Patient, timepoints pooled | `~ condition` |
| Timepoint within a condition | Within a patient | Sample, patient blocks | `~ patient_id + timepoint` |
| Case vs control at one timepoint | Between patients | Sample, one per patient | `~ condition` |
| Condition by timepoint | Both | Sample, patient random | `~ condition * timepoint + (1 \| patient_id)` |

## The four families

**Family 1, condition.** The headline contrast. Condition varies between
patients, so the patient is the unit and all their timepoints pool into one
profile.

**Family 2, within-condition timepoints, paired.** Every pairwise timepoint
contrast inside each condition. Only patients sampled at both timepoints of a
pair are kept, which is what makes `~ patient_id + timepoint` a genuine paired
design. A patient present at one timepoint only cannot inform the difference,
and leaving them in makes the design matrix rank deficient in exactly the cases
that matter.

The script re-checks completeness after the minimum-cell filter, because a
patient can lose one of their two samples to that filter and break the pairing
without anything else noticing.

**Family 3, between-condition at each timepoint, unpaired.** Each patient
contributes one sample per timepoint, so there is nothing to block on.
Timepoints where either arm has fewer than `MIN_UNITS_PER_GROUP` patients are
skipped by a guard rather than fitted.

**Family 4, the interaction.** The first three never ask whether the case arm
changes differently from control. Family 1 collapses time, family 2 looks
inside one condition, family 3 freezes a timepoint. The interaction is the
difference of those differences, which is the trajectory question.

## Why the interaction needs a different tool

DESeq2 cannot fit it. Patient is nested inside condition, so
`~ condition * timepoint + patient_id` has a rank-deficient design matrix: a
patient's fixed effect already encodes their condition.

Patient has to enter as a random intercept, which means a mixed model per gene.
[`variancePartition::dream`](https://doi.org/10.1093/bioinformatics/btaa687)
does exactly this: `voomWithDreamWeights` computes precision weights under the
mixed model, `dream` fits the per-gene model, `eBayes` moderates variances
across genes the way limma does.

Two outputs:

- **Omnibus**, a moderated F over all interaction coefficients jointly. Does
  the trajectory differ anywhere. Written to its own subfolder, because an F
  statistic is unsigned and cannot rank a GSEA list.
- **Per-timepoint**, one signed coefficient each. How much more the arms differ
  at that timepoint than they already did at the reference.

## The guards

Every one of these returns `NULL` rather than a bad table, and the skip is
logged:

```r
if (ncol(pb) < 2) return(NULL)                        # not enough units
if (min(n_test, n_ref) < MIN_UNITS_PER_GROUP) NULL    # no within-group dispersion
if (qr(mm)$rank < ncol(mm)) return(NULL)              # singular or nested design
if (nrow(mm) - ncol(mm) < 1) return(NULL)             # no residual df
if (sum(keep) < MIN_GENES_TESTED) return(NULL)        # fit not worth reporting
```

A pseudobulk unit with fewer than `MIN_CELLS` cells is dropped before any of
this. Summing 3 cells produces a profile that is mostly noise, and it would
enter the model with the same weight as one built from 3,000.

## The run manifest

Every attempted contrast gets a row, with the number of units on each side and
a status:

| Status | Meaning |
|:--|:--|
| `ok` | Fitted, at least 3 units per arm |
| `exploratory` | Fitted, fewer than 3 units in one arm |
| `skipped` | A guard fired, with the reason |

This is the sheet to point reviewers at. A comparison that was skipped is
distinguishable from one that was never asked for, and a result resting on two
units per arm is labelled rather than read as though it had the same standing
as the headline contrast.

## Two statistics, two jobs

| Column | What it is | Used for |
|:--|:--|:--|
| `avg_log2FC` | `lfcShrink` with `ashr` | Volcano plots. Effect size. |
| `stat` | Unshrunken Wald statistic | GSEA ranking. Evidence. |

Both are written to every table. Shrinkage pulls noisy low-count estimates
toward zero, which is right for a plot about effect size and wrong for a
ranking, because it discards the information about how well-estimated each
effect was.

## Biotype annotation

DEG tables are annotated in place against a GENCODE GTF, flagging pseudogenes,
lncRNAs, immunoglobulin and TCR genes, and uncharacterised entries with
`remove_flag`.

Annotating in place rather than writing a second file means GSEA reads the same
csv, and there is one table per contrast rather than two that can drift. The
flag is not applied automatically: volcanoes filter on it, and GSEA has
`GSEA_DROP_FLAGGED` off by default.
