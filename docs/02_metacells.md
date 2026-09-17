# What a metacell is

## The problem

A droplet gives a few thousand UMIs from a cell holding a few hundred thousand
mRNA molecules. Most genes in most cells read zero, and a zero does not
distinguish "off" from "expressed, not sampled".

What the data does have is redundancy. In 250,000 PBMCs there are not 250,000
distinct transcriptional states. There are perhaps a few thousand, each
represented by tens to hundreds of cells that differ from each other by
sampling noise. Identify those groups, sum their counts, and the result is one
deep profile per state instead of many shallow ones.

That is a metacell: a group of cells whose differences are consistent with
being independent draws from the same distribution.

## Metacell against cluster

Confusing the two causes real mistakes.

**A cluster** is meant to be a cell type or state. There are 10 to 30. Each
holds cells that genuinely differ from each other, and how much difference is
acceptable is a resolution parameter turned by hand.

**A metacell** is meant to be a sampling unit. There are thousands. Each holds
cells that are not meaningfully different, by construction, and the algorithm's
job is to stop growing a group before it absorbs real variation. Size is set by
the UMI depth needed, not by the number of expected cell types.

Metacells are therefore upstream of clustering. Here they are grouped into
blocks afterwards ([07](07_annotation.md)), and the blocks carry the labels.

## How MC2 builds them

[Metacell-2](https://doi.org/10.1186/s13059-022-02667-1), in five steps.

**1. Feature selection.** Keep genes whose variance across cells exceeds what
Poisson sampling alone would produce. Only these decide which cells are
similar. Everything else is still counted, it just gets no vote.

**2. Similarity graph.** A balanced k-nearest-neighbour graph over cells on
those features. Balanced means an edge survives only if both cells consider the
other a close neighbour, which stops one high-count cell becoming everybody's
neighbour.

**3. Partition.** Cut the graph into groups sized to hit a target cell count
and a target UMI count. Cells that fit no coherent group are marked outliers
rather than forced into the nearest one.

**4. Divide and conquer.** All of the above at once on 250,000 cells is not
tractable, so the data is split into piles of a few thousand cells, metacells
are computed within each pile, and the results are reconciled across piles.
Pile boundaries are a computational device; the reconciliation step is what
keeps them out of the output.

**5. Collect.** Sum raw counts within each metacell. The output object has one
row per metacell, at a UMI depth in the hundreds of thousands.

## The costs

**Cells cannot be counted by counting metacells.** MC2 sizes metacells to a
target cell count, so the number of metacells in a population is not
proportional to that population's abundance. Every composition analysis here
counts cells, each carrying its metacell's label. Counting metacells would be
wrong, and wrong invisibly. See [09](09_composition.md).

**Rare populations can vanish.** A state with 30 cells across the dataset will
not form its own metacell at a target size of 72. It gets absorbed, or called
outlier. MC2's rare gene module detector exists for this, and steps 02 and 04
both run it and report what it found.

**Outliers are real data.** Typically 5 to 10 percent of cells end up
unassigned. A rising outlier rate on a subset usually means
`target_metacell_size` is too large for how much structure is left, which makes
it a diagnostic rather than a nuisance.

## When it is the wrong tool

Metacells help when there are many cells, shallow coverage, and a question
about gene-level detail. They do not help with few cells, with questions about
individual cells, or with populations that are genuinely rare rather than
merely small.

For differential expression, pseudobulk is required regardless, because the
replicate is the donor and not the cell ([10](10_differential_expression.md)).
