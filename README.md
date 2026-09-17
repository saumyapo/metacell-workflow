# Metacell Workflow (MC2)

![R](https://img.shields.io/badge/R-4.4-276DC3?style=flat-square&logo=r&logoColor=white)
![Python](https://img.shields.io/badge/Python-3.10-3776AB?style=flat-square&logo=python&logoColor=white)
![Seurat](https://img.shields.io/badge/Seurat-v5-black?style=flat-square)
![MC2](https://img.shields.io/badge/metacells-0.9.5-8B1A1A?style=flat-square)
[![license](https://img.shields.io/badge/license-MIT-green?style=flat-square)](LICENSE)

End-to-end metacell pipeline for single-cell RNA-seq with paired V(D)J: raw
matrices in, cell type annotations, composition statistics, differential
expression, pathway enrichment and clonal repertoire out.

## What a metacell is, in one paragraph

A droplet gives a few thousand UMIs from a cell holding a few hundred thousand
mRNA molecules, so most genes in most cells read zero and a zero is ambiguous.
Metacells ([Ben-Kiki et al. 2022](https://doi.org/10.1186/s13059-022-02667-1))
pool cells that are statistically indistinguishable into single profiles of a
few hundred thousand UMIs each. The result is deep enough to see genes that
were dropout noise at the single-cell level, and still fine grained enough that
real states are not averaged away.

The pipeline does that twice: once across the whole dataset to get broad
lineages, then again inside one lineage to resolve states. Composition,
differential expression, GSEA and clonal repertoire all hang off those two
annotations.

New to any of this, start with [docs/01_background.md](docs/01_background.md).

## The pipeline

```
   raw 10x h5              +------------------------------------------+
        |                  |  01  preprocess                          |
        v                  |      dedupe genes, TCR doublet filter,   |
   merged .h5ad  <---------+      attach design, per-sample QC        |
        |                  +------------------------------------------+
        v
   +--------------------+          +----------------------------+
   | 02  global         |          | 03  global annotation      |
   |     metacells      |--------->|     gene modules, block    |
   |     MC2 divide     |          |     heatmap, lineage calls |
   |     and conquer    |          +-----------+----------------+
   +--------------------+                      |  broad_label
                                               v
   +--------------------+          +----------------------------+
   | 04  subset         |          | 05  subset annotation      |
   |     metacells      |--------->|     finer modules, state   |
   |     purity gate,   |          |     calls                  |
   |     finer grain    |          +-----------+----------------+
   +--------------------+                      |  cluster_label
                                               v
        +--------------------------------------+-------------------------+
        v                    v                  v                        v
  06 stacked bars     08 pseudobulk DE    10 clonal repertoire   12 external labels
  07 CLR models       09 GSEA             11 clonal statistics
```

| Step | File | Lang | What it does |
|:--|:--|:--|:--|
| 01 | `scripts/01_preprocess.py` | Py | Load per-sample matrices, collapse duplicated gene symbols, remove dual-TRB doublets, attach the design, merge, write QC |
| 02 | `scripts/02_global_metacells.ipynb` | Py | Cell and gene QC, excluded/lateral/noisy gene lists, MC2 divide and conquer over the whole dataset |
| 03 | `scripts/03_global_modules_annotation.Rmd` | R | Gene-gene correlation modules, module-score heatmap, metacell blocks, lineage calls |
| 04 | `scripts/04_subset_metacells.ipynb` | Py | Pull one lineage out, purity gate, strip stale MC2 fields, rerun MC2 at finer grain |
| 05 | `scripts/05_subset_modules_annotation.Rmd` | R | The same module machinery inside one lineage, state calls |
| 06 | `scripts/06_composition_stacked_bars.Rmd` | R | Per-sample composition, donor-weighted stacked and dodged bars |
| 07 | `scripts/07_composition_clr_models.Rmd` | R | CLR transform, mixed models across timepoints, forest plots |
| 08 | `scripts/08_pseudobulk_degs.Rmd` | R | Four DESeq2 contrast families plus a `dream` mixed-model interaction |
| 09 | `scripts/09_pseudobulk_gsea.Rmd` | R | GSEA ranked on the Wald statistic, one xlsx per DEG csv |
| 10 | `scripts/10_clonal_repertoire.Rmd` | R | scRepertoire, clone size, occurrence spectra, clone landscape |
| 11 | `scripts/11_clonal_statistics.Rmd` | R | Rarefied diversity, Renyi profiles, Gini, Tversky sharing, Startrac, tests |
| 12 | `scripts/12_external_label_mapping.Rmd` | R | Transfer an external per-cell binary call onto the annotated object |

Steps 04 and 05 are rerun once per lineage. Everything from 06 onwards runs on
either the global or a subset object, selected by one `LEVEL` variable at the
top of each file.

## Quick start

```bash
git clone https://github.com/saumyapo/metacell-workflow.git && cd metacell-workflow
mamba env create -f env/metacell.yml && conda activate metacell   # or use the container
$EDITOR config/samples.csv config/vdj_library_map.csv config/config.yaml
python3 scripts/01_preprocess.py
jupyter lab scripts/02_global_metacells.ipynb
```

Then work down the table. Every step reads what the previous one wrote, and
nothing outside `config/` needs editing between studies.

## Configuration

No sample IDs and no absolute paths anywhere in `scripts/`. Four files hold the
study description and every path:

| File | Holds |
|:--|:--|
| `config/samples.csv` | One row per sample: patient, timepoint, condition, plotting label |
| `config/vdj_library_map.csv` | Library folder token to sample ID, because library naming is rarely derivable by pattern |
| `config/config.yaml` | Paths and every tunable for the Python steps |
| `config/config.R` | The same paths and factor levels for the R steps |

Pointing the pipeline at a different study is two CSVs and a path. Tuning knobs
that are read while looking at output (`NCLUST`, `column_km`, the HVG
thresholds) are declared at the top of the script that uses them, since they are
set from a diagnostic printed in that same file. See
[docs/06_parameters.md](docs/06_parameters.md).

## Environments

The Python and R stacks are separate images, because their system library
requirements conflict and an R package bump should not invalidate the MC2
image.

```bash
apptainer build env/metacell.sif   env/metacell.def
apptainer build env/metacell-r.sif env/metacell-r.def

apptainer exec --bind "$PWD":"$PWD" --pwd "$PWD" env/metacell.sif python3 scripts/01_preprocess.py
```

The Python image is generated from the same `env/metacell.yml` a local conda
env would use, so the container and the bare env cannot drift apart. Building
needs root or fakeroot, so a login node or a CI runner rather than a compute
node. Details in [docs/14_reproducibility.md](docs/14_reproducibility.md).

## Documentation

[`docs/`](docs/) is the reasoning behind the code: what a metacell is, how to
set `target_metacell_size`, why the lateral gene list matters more than it
looks, how to read a block heatmap, and what every printed diagnostic means.
Start at [docs/README.md](docs/README.md).

## Citing

> Ben-Kiki O, Bercovich A, Lifshitz A, Tanay A. Metacell-2: a divide-and-conquer
> metacell algorithm for scalable scRNA-seq analysis. *Genome Biology* 23, 100 (2022).
