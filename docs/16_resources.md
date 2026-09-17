# Resources

External reading. The pipeline docs explain the choices made here; these
explain the underlying methods properly, which is worth doing once.

## Starting from scratch in scRNA-seq

- [Single-cell best practices](https://www.sc-best-practices.org/) (Heumos et
  al.). The current reference text. Chapters on QC, normalisation, annotation
  and compositional analysis are all directly relevant here.
- [Orchestrating Single-Cell Analysis with Bioconductor](https://bioconductor.org/books/release/OSCA/)
  (OSCA). Same ground from the Bioconductor side, with more statistical detail.
- [Current best practices in single-cell RNA-seq analysis: a tutorial](https://doi.org/10.15252/msb.20188746)
  (Luecken and Theis 2019). Older, still the clearest single overview of what
  the steps are and why.

## Metacells

- [Metacell-2 paper](https://doi.org/10.1186/s13059-022-02667-1) (Ben-Kiki et
  al. 2022). The algorithm this pipeline runs.
- [metacells package docs](https://metacells.readthedocs.io/en/latest/) and
  [API reference](https://metacells.readthedocs.io/en/latest/API.html).
- [Official vignettes](https://github.com/tanaylab/metacells-vignettes). Worth
  running once end to end on their data before running this pipeline on new
  data.
- [MetaCell (original R implementation)](https://doi.org/10.1186/s13059-019-1812-2)
  (Baran et al. 2019). Background only; superseded by MC2.
- [SuperCell](https://doi.org/10.1186/s12859-022-04861-1) (Bilous et al. 2022).
  An alternative implementation, and the benchmark the graining level in
  [06](06_parameters.md) is set from.
- [SEACells](https://doi.org/10.1038/s41587-023-01716-9) (Persad et al. 2023).
  Another alternative, archetype-based, works on ATAC too.

## QC and doublets

- [miQC](https://doi.org/10.1371/journal.pcbi.1009290) (Hippen et al. 2021).
  Mixture-model mitochondrial filtering, reimplemented and left commented out
  in steps 02 and 04.
- [Benchmarking computational doublet-detection methods](https://doi.org/10.1016/j.cels.2020.11.008)
  (Xi and Li 2021).
- [scirpy](https://doi.org/10.1093/bioinformatics/btaa611) (Sturm et al. 2020).
  The chain-QC logic behind the dual-TRB doublet filter in step 01.

## Annotation

- [ComplexHeatmap reference book](https://jokergoo.github.io/ComplexHeatmap-reference/book/).
  The block heatmaps in 03 and 05 are built on it, and `column_km` is the
  resolution knob.
- [Azimuth PBMC references](https://azimuth.hubmapconsortium.org/). Useful as a
  marker-panel sanity check, including for deciding which markers cannot
  separate two states and therefore belong in `panel` rather than `sig`.

## Composition

- [Compositional analysis chapter, sc-best-practices](https://www.sc-best-practices.org/conditions/compositional.html).
  The clearest explanation of why raw proportions cannot be tested directly.
- [propeller](https://doi.org/10.1093/bioinformatics/btac582) (Phipson et al.
  2022) and [scCODA](https://doi.org/10.1038/s41467-021-27150-6) (Büttner et
  al. 2021). Two heavier alternatives to the CLR plus mixed model used here.

## Differential expression

- [Confronting false discoveries in single-cell differential expression](https://doi.org/10.1038/s41467-021-25960-2)
  (Squair et al. 2021). The paper to read if pseudobulk still looks like an
  optional step.
- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).
  Design matrices, contrasts and shrinkage.
- [dream](https://doi.org/10.1093/bioinformatics/btaa687) (Hoffman and Roussos
  2021). The mixed-model route used for the interaction family.

## Pathway enrichment

- [fgsea vignette](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html).
- [clusterProfiler book](https://yulab-smu.top/biomedical-knowledge-mining-book/).
- [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) for what is actually in
  Hallmark, GO BP, KEGG and Reactome.

## Clonal repertoire

- [scRepertoire documentation](https://www.borch.dev/uploads/screpertoire/).
  Vignettes for every function used in steps 10 and 11.
- [rnabio.org clonality module](https://rnabio.org/module-08-scrna/0008/06/01/Clonality/).
  The workflow step 10 is structured after.
- [Mhanna et al. 2024](https://doi.org/10.1016/j.crmeth.2024.100738).
  Occurrence-based clone size bins, Renyi profiles, rarefaction curves.
- [10x V(D)J algorithm documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-5p-vdj-algorithm).
  What `filtered_contig_annotations.csv` columns actually mean.

## Tools

- [Seurat v5](https://satijalab.org/seurat/)
- [scanpy](https://scanpy.readthedocs.io/) and
  [anndata](https://anndata.readthedocs.io/)
- [Apptainer](https://apptainer.org/docs/user/latest/) for the containers in
  [14](14_reproducibility.md)
