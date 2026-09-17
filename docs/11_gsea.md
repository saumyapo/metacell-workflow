# Pathway enrichment

Step 09.

## Ranking statistic

Genes are ranked by the DESeq2 Wald statistic, or by the `dream` moderated t
for the interaction family. Not by fold change.

The Wald statistic is the effect divided by its standard error. A large but
noisy fold change in a low-count gene does not outrank a moderate,
well-estimated one, which is what a ranked enrichment test needs.

Ranking on log2FC alone puts low-expression genes at both ends of the list,
because that is where the noisiest estimates are, and the enrichment follows
them there. The shrunken `avg_log2FC` in the same tables is for volcanoes,
where the question is effect size rather than evidence.

`RANK_COL` falls back to `avg_log2FC` if `stat` is absent, so a table from a
different tool still runs, but the fallback is recorded in the manifest.

## Duplicate symbols

`fgsea` indexes the ranked vector by name, so a duplicated symbol silently
overwrites its twin and the set membership lookup goes wrong.

```r
d %>% arrange(desc(abs(.data[[rank_col]]))) %>% distinct(gene, .keep_all = TRUE)
```

Duplicates are resolved to the strongest statistic before ranking rather than
left to whichever row happened to come last.

## Gene sets

Hallmark, GO Biological Process, KEGG legacy and Reactome, from `msigdbr`.

`minGSSize = 15`, `maxGSSize = 500`. The lower bound keeps out sets too small
to have a stable enrichment score; the upper bound keeps out sets so broad they
are enriched for almost any list.

`pvalueCutoff = 1` so nothing is filtered at run time. Filtering happens in a
separate pass that reads the finished workbooks back, which means the
unfiltered result is always on disk and the significance threshold can change
without rerunning the enrichment.

## The omnibus family

Omnibus interaction tables carry `F_stat` and no signed statistic, so they
cannot be ranked. They are scanned anyway and logged as `not rankable`, so
every DEG csv is accounted for in the manifest rather than silently absent.

That distinction, "looked and it cannot be done" against "not there", is the
point of the manifest.

## Filename grammar

```
FILE_PREFIX__unit__contrast.csv          ->  FILE_PREFIX__unit__contrast__GSEA.xlsx
```

One DEG csv, one GSEA xlsx, same stem, matching family folder. The manifest is
a complete one-to-one map between them, and a mismatch in row counts is the
first thing to check when something is missing.

## Filtering

A second pass reads the finished workbooks and writes `__GSEA_filtered.xlsx`
alongside each, keeping `|NES| > 1.5` and `q < 0.05`.

An empty filtered file is written when nothing passes, so every unfiltered file
keeps a one-to-one partner and a missing file always means a genuine failure
rather than a null result.

> [!TIP]
> `read.xlsx` returns `NULL` for a header-only sheet and drops all-NA columns
> by default, both of which break `bind_rows` downstream. The `read_gsea()`
> helper normalises this once so the rest of the code does not have to guard
> against it.
