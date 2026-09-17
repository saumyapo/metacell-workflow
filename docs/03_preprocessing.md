# Preprocessing

Step 01 (`scripts/01_preprocess.py`). Nothing here is metacell-specific. It
produces the merged object that step 02 picks up.

## Duplicated gene symbols

Cell Ranger references contain genes that map to the same symbol more than
once. Left alone they stay separate columns, so a cell's expression of that
gene is split across them and every downstream sum is wrong by an unpredictable
amount.

Step 01 sums the duplicated columns into one and keeps the first copy's
metadata. This has to happen before anything else, because every filter and
total downstream depends on it. The check afterwards prints an empty vector per
sample when it worked.

## Doublets

Two mechanisms, catching different things.

### Expression-based (Scrublet)

Simulates artificial doublets by adding pairs of real profiles, then scores
each barcode on how much it looks like one. Available behind `RUN_SCRUBLET`,
off by default.

Its weakness is structural: a doublet of two cells of the same type looks like
one cell of that type with more UMIs. Homotypic doublets are invisible to it,
and in a T-cell-dominated PBMC sample most doublets are homotypic.

### V(D)J-based (dual productive TRB)

This is what step 01 runs by default.

A 10x barcode labels a GEM, not a cell, and the gene expression and V(D)J
libraries are barcoded from the same partition
([Zheng et al. 2017](https://doi.org/10.1038/ncomms14049)). TCR beta allelic
exclusion is near-complete, so a single T cell expresses one productive beta
chain ([Brady et al. 2010](https://doi.org/10.4049/jimmunol.1001158)). Two
distinct productive betas on one barcode therefore mean two T cells in one
droplet.

That is a physical observation about the droplet rather than an inference from
expression, which is why it catches the homotypic T-T doublets Scrublet cannot.
It is the same principle as the multichain category in
[scirpy's chain QC](https://doi.org/10.1093/bioinformatics/btaa611).

**Alpha is deliberately not used.** TCR alpha allelic exclusion is leaky, and 10
to 30 percent of genuine T cells carry two productive alpha chains
([Padovan et al. 1993](https://doi.org/10.1126/science.8493531); reviewed in
[Schuldt and Binstadt 2019](https://doi.org/10.4049/jimmunol.1801430)). Gating
on dual alpha would discard real cells at scale. The script counts them,
reports them, and keeps them.

**Ambient TCR is filtered out.** mRNA from lysed T cells gets assembled into
low-support contigs on unrelated barcodes
([Young and Behjati 2020](https://doi.org/10.1093/gigascience/giaa151);
[Fleming et al. 2023](https://doi.org/10.1038/s41592-023-01943-7)). So the
second beta has to clear both an absolute UMI floor and a fraction of the
dominant beta:

```yaml
min_second_beta_umis: 2      # absolute floor
min_second_beta_frac: 0.30   # and at least 30% of the top-ranked TRB
```

Both live in `config/config.yaml`. Raising either makes the filter more
conservative and catches fewer doublets.

### Scope, and what belongs in methods

This removes only doublets where at least one partner is a T cell and a TCR was
recovered. It complements expression-based calling rather than replacing it,
and current benchmarking recommends keeping both
([Xi and Li 2021](https://doi.org/10.1016/j.cels.2020.11.008);
[Heumos et al. 2023](https://doi.org/10.1038/s41576-023-00586-w)).

The script prints a sanity check against the 10x expectation of roughly 0.8
percent multiplets per 1,000 recovered cells. The observed TCR-doublet rate
should sit well below that, since only T-containing doublets with a recovered
TCR are visible.

> [!WARNING]
> If composition claims rest on exact cell counts, turn Scrublet on. Two
> filters change the denominator of every frequency in the paper. Decide once,
> early, and record which was used.

## Design metadata

The whole design comes from `config/samples.csv` via a join on `sample_origin`.
Nothing is parsed out of the sample ID string.

Parsing IDs with a regex works until one sample is named slightly differently,
at which point it fails silently or matches the wrong pattern. A join either
finds the row or errors.

The V(D)J library map is a separate CSV for the same reason: library folder
names are rarely derivable from sample IDs by pattern, particularly when visit
numbering is inconsistent across donors.

## QC output

One workbook, appended to by step 02:

| Sheet | Contents |
|:--|:--|
| `preprocess_qc_summary` | Cells and genes loaded, UMI stats, TCR recovery, doublets removed, retention, per sample plus a total row |
| `umi_percentiles` | Full UMI distribution per sample |
| `excluded_per_sample` | Written by 02: attrition through each filter, plus outlier rate |

Read the retention column across samples before anything else. If one sample
loses far more cells than the others, that is a batch problem, and every
proportion downstream inherits it.
