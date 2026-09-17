# Background

Start here if metacells are new. This page covers what problem they solve,
where they sit relative to the more familiar steps of an scRNA-seq analysis,
and the vocabulary the rest of the docs uses without explaining.

## The data

A droplet scRNA-seq run (10x Chromium and similar) partitions cells into
nanolitre droplets, barcodes the mRNA inside each one, and sequences it. After
Cell Ranger, the result is a counts matrix: cells on one axis, genes on the
other, and a UMI count in each entry. One UMI is one captured mRNA molecule.

Two numbers set everything that follows. A typical cell contains a few hundred
thousand mRNA molecules. A typical droplet yields a few thousand UMIs. So each
cell is sampled at roughly one percent, and the sampling is multinomial from a
pool where a few dozen genes take most of the probability mass.

The consequence is that most entries in the matrix are zero, and a zero is
ambiguous. It can mean the gene is off, or it can mean the gene is expressed at
a moderate level and was not drawn. Normalisation does not fix this. The
information is not in the matrix to begin with.

## Where coarse-graining comes in

The usual response is to group cells and work with the group. Three versions of
this appear in most pipelines, and they are not interchangeable:

**Clustering.** Groups cells into 10 to 30 units meant to be cell types or
states. The groups are large and internally heterogeneous on purpose, and how
heterogeneous is a resolution parameter set by hand.

**Pseudobulk.** Sums counts within a sample and a cell type to get one profile
per biological replicate. This is a statistical necessity for differential
expression (see [10](10_differential_expression.md)), not a way to see more
structure. It collapses all within-sample variation.

**Metacells.** Groups cells into thousands of small units that are meant to be
internally homogeneous, so that summing their counts gives a deep profile of
one transcriptional state rather than a blend of several. Size is set by how
many UMIs are needed, not by how many cell types are expected.

Metacells sit upstream of clustering, not in place of it. In this pipeline the
metacells are grouped into blocks afterwards, and the blocks get the cell type
labels. Pseudobulk still happens later, for the same statistical reasons it
always did.

## Which implementation

Several tools build metacell-like units. This pipeline uses Metacell-2 (MC2),
the Python package from the Tanay lab
([Ben-Kiki et al. 2022](https://doi.org/10.1186/s13059-022-02667-1)).

| Tool | Notes |
|:--|:--|
| [Metacell-2 (MC2)](https://github.com/tanaylab/metacells) | Python. Divide and conquer, so it scales to millions of cells. Explicit outlier handling, explicit lateral gene support. Used here. |
| [MetaCell (R)](https://github.com/tanaylab/metacell) | The original. Superseded by MC2 and not maintained for new work. |
| [SuperCell](https://doi.org/10.1186/s12859-022-04861-1) | R, walktrap on a kNN graph. Simpler and faster. The same paper is the benchmark that this pipeline's graining level is set from. |
| [SEACells](https://doi.org/10.1038/s41587-023-01716-9) | Python, archetypal analysis. Works on ATAC as well as RNA. |

MC2 is the choice here mainly for the outlier and lateral gene handling, which
are what make the two-pass structure in [08](08_two_passes.md) work.

## Where metacells do not help

- Few cells. Pooling needs redundancy to pool.
- Questions about individual cells: doublet calling, trajectory inference at
  single-cell resolution.
- Genuinely rare populations. A state with 30 cells in the whole dataset will
  not form its own metacell. MC2 has a separate rare gene module detector for
  this, covered in [04](04_building_metacells.md).
- Differential expression on its own. Pseudobulk is still required, because the
  replicate is the donor. Metacells buy a cleaner annotation to pseudobulk
  within, not a way around the statistics.

## Vocabulary

| Term | Meaning |
|:--|:--|
| Metacell | A group of cells whose differences are consistent with sampling noise alone. Counts are summed within it. |
| Graining level (gamma) | Cells divided by metacells. Set by `target_metacell_size`. See [06](06_parameters.md). |
| Pile | A chunk of a few thousand cells that MC2 processes independently during divide and conquer. Computational, not biological. |
| Outlier | A cell MC2 could not place in any metacell. Gets `metacell = -1`. Normally 5 to 10 percent. |
| Excluded gene | Dropped from the object entirely. See [05](05_gene_lists.md). |
| Lateral gene | Counted, but not allowed to influence cell-cell similarity. |
| Noisy gene | Given extra slack in outlier detection because it is bursty. |
| Gene module | A set of genes that co-vary across metacells, found by correlation clustering in [07](07_annotation.md). |
| Block | A k-means group of metacells on the module score matrix. Blocks carry the cell type labels. |
| Pseudobulk | Counts summed to one profile per replicate unit before testing. |
| CLR | Centred log-ratio. Takes proportions out of the simplex so linear models apply. See [09](09_composition.md). |
| Clone | A set of T cells carrying the same receptor rearrangement. See [12](12_clonal.md). |

## What to read next

[02 What a metacell is](02_metacells.md) for the algorithm. Papers, tutorials
and package documentation are collected in [16 Resources](16_resources.md).
