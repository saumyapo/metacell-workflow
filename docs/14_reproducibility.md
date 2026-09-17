# Reproducibility

## Two environments, on purpose

The Python side (steps 01, 02, 04) and the R side (03, 05 to 12) are separate
images.

They have different, occasionally conflicting, system library requirements, and
splitting them means an R package bump cannot invalidate the MC2 image. The one
place they meet is `reticulate` in steps 03 and 05, which reads the `.h5ad`
files. That is handled by binding the Python environment into the R container
rather than installing conda twice.

## Python

`env/metacell.yml` is a full conda export with build strings, which is exact
and Linux-x86-64 only. The container is built from the same file, so a local
env and the image cannot drift.

```bash
mamba env create -f env/metacell.yml && conda activate metacell
```

Core versions: `metacells 0.9.5`, `scanpy 1.10.4`, `anndata 0.10.9`,
`numpy 2.2.x`, Python 3.10.

> [!WARNING]
> The export pins `numpy` twice: `numpy-base=2.2.5` from conda and
> `numpy==2.2.6` from pip. Pip wins at install time and the result works, but a
> pip package overwriting a conda one breaks on some later solve. Resolve it to
> a single source before building an image intended for citation.

## R

`env/metacell-r.def` builds on `rocker/r-ver:4.4.2` with a Posit Package
Manager snapshot date pinned in `Rprofile.site`:

```r
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2025-01-15"))
```

One line decides what "current" means for every CRAN package, which is easier
to reason about, and to bump deliberately, than a list of pinned versions.
Bioconductor is pinned by release (3.20), which pins its own set.

Two hard requirements the scripts assert on:

```r
stopifnot(packageVersion("Seurat") >= "5.0.0")
stopifnot(packageVersion("scRepertoire") >= "2.0.0")
```

scRepertoire 1.x and 2.x differ enough that the clonal notebooks will not run
on 1.x, and failing at the top beats failing in the middle.

## Containers

Build once, from the definition files:

```bash
apptainer build env/metacell.sif   env/metacell.def
apptainer build env/metacell-r.sif env/metacell-r.def
```

Building needs root or fakeroot, which most clusters do not allow on compute
nodes. A login node or a CI runner works. If the images are pushed to a
registry, `apptainer pull oras://<registry>/<image>:<tag>` replaces the build
step on the cluster entirely.

Running them:

```bash
# python steps
apptainer exec --bind "$PWD":"$PWD" --pwd "$PWD" env/metacell.sif python3 scripts/01_preprocess.py

# R steps, with the python env bound in for reticulate
apptainer exec --bind "$PWD":"$PWD" --bind /path/to/envs/metacell:/opt/pyenv --pwd "$PWD" --env RETICULATE_PYTHON=/opt/pyenv/bin/python env/metacell-r.sif Rscript -e 'rmarkdown::render("scripts/03_global_modules_annotation.Rmd")'
```

`--bind "$PWD":"$PWD" --pwd "$PWD"` keeps paths identical inside and outside the
container, which is what lets `config/config.R` and `config.yaml` work unchanged
either way.

Both images record what they contain at build time: `/opt/metacell.lock` and
`/opt/metacell.pip.txt` for the Python image, `/opt/r-packages.csv` for the R
image. A container can therefore be audited without running it.

## Seeds

`random_seed: 123` goes to every MC2 call, `set.seed(42)` before the heatmap
clustering in 03 and 05, and GSEA uses `seed = TRUE` with `set.seed(2024)`.

This makes a rerun on the same input reproduce. It does not make results
robust: if a conclusion changes when the seed changes, the seed is not what
fixed it.

## What is not pinned

Cell Ranger version and reference, GENCODE annotation for the biotype step, and
MSigDB collection version (whatever `msigdbr` ships). These belong in the
methods section; the pipeline does not capture them.

## Session info

Every Rmd prints package versions in its header, and the knitted HTML carries
them. Keeping the HTML alongside the outputs is the cheapest provenance record
available: it says what ran, with what, and what it printed.
