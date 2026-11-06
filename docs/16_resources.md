<div align="right">
  <small><em>Author: Saumya Pothukuchi & Claude AI - Opus 4.8 (Anthropic, 2026)</em></small>
</div>

# Resources

* The pages in this folder explain the choices made in this pipeline. The resources below explain the underlying methods in full, and are worth reading once for any step that is new.
* Each entry says what it is and where it connects to this pipeline. Entries are grouped in pipeline order.

## Where to start

| If you are | Read first |
|:--|:--|
| New to scRNA-seq analysis | Single-cell best practices, chapters on QC, annotation and conditions |
| Comfortable with scRNA-seq, new to metacells | The MC2 paper, then the Bilous et al. 2024 review, then [02](02_metacells.md) |
| Deciding between the rerun and reuse subset paths | The MC2 paper (divide and conquer section), the official vignettes, then [08b](08b_rerun_or_reuse.md) |
| Reviewing the statistics | Squair et al. 2021, the compositional chapter of sc-best-practices, the `dream` paper |

## General scRNA-seq

- [Single-cell best practices](https://www.sc-best-practices.org/) (Heumos et al.). The current reference text, maintained online. The QC, annotation, differential gene expression and compositional analysis chapters all map directly onto steps here.
- [Best practices for single-cell analysis across modalities](https://doi.org/10.1038/s41576-023-00586-w) (Heumos et al. 2023). The peer-reviewed summary of the book, with recommendations per step.
- [Orchestrating Single-Cell Analysis with Bioconductor](https://bioconductor.org/books/release/OSCA/) (OSCA). The same ground from the Bioconductor side, with more statistical detail.
- [Current best practices in single-cell RNA-seq analysis: a tutorial](https://doi.org/10.15252/msb.20188746) (Luecken and Theis 2019). Older, still the clearest short overview of what each step does and why.

## Metacells

- [Metacell-2: a divide-and-conquer metacell algorithm](https://doi.org/10.1186/s13059-022-02667-1) (Ben-Kiki et al. 2022). The algorithm this pipeline runs. The divide and conquer section explains how piles are regrouped by similarity in the final phase, which is central to the rerun or reuse question ([08b](08b_rerun_or_reuse.md)).
- [metacells documentation](https://metacells.readthedocs.io/en/latest/) and [API reference](https://metacells.readthedocs.io/en/latest/API.html). Every function called in steps 02 and 04, with its parameters and defaults.
- [Official vignettes](https://github.com/tanaylab/metacells-vignettes) (Tanay lab). One-pass, iterative and projection workflows. The iterative vignette shows how lateral and excluded gene lists are refined across repeated runs on the whole dataset. Worth running once on the example data before running this pipeline on new data.
- [MCView](https://github.com/tanaylab/MCView) (Tanay lab). An interactive viewer for MC2 output, useful for exploring metacells and gene-gene relationships alongside the static heatmaps here.
- [Building and analyzing metacells in single-cell genomics data](https://doi.org/10.1038/s44320-024-00045-6) (Bilous et al. 2024). A review and tutorial covering when metacells help, when they do not, quality metrics, and analysis at the metacell level. Accompanied by a [tutorial](https://github.com/GfellerLab/MetacellAnalysisTutorial) and a [toolkit](https://github.com/GfellerLab/MetacellAnalysisToolkit).
- [Metacells untangle large and complex single-cell transcriptome networks](https://doi.org/10.1186/s12859-022-04861-1) (Bilous et al. 2022). SuperCell, an alternative implementation, and the benchmark of graining levels that [06](06_parameters.md) uses.
- [SEACells](https://doi.org/10.1038/s41587-023-01716-9) (Persad et al. 2023). Archetype-based metacells, also applicable to ATAC-seq.
- [MetaCell, the original R implementation](https://doi.org/10.1186/s13059-019-1812-2) (Baran et al. 2019). Background only; superseded by MC2.

## Preprocessing, QC and doublets

- [Massively parallel digital transcriptional profiling of single cells](https://doi.org/10.1038/ncomms14049) (Zheng et al. 2017). The 10x barcoding design: why one barcode marks a droplet, not a cell.
- [scirpy](https://doi.org/10.1093/bioinformatics/btaa611) (Sturm et al. 2020). Receptor chain QC, including the multichain category behind the dual-TRB doublet filter in step 01.
- [Benchmarking computational doublet-detection methods](https://doi.org/10.1016/j.cels.2020.11.008) (Xi and Li 2021). How expression-based doublet callers compare, and what they miss.
- [miQC](https://doi.org/10.1371/journal.pcbi.1009290) (Hippen et al. 2021). Mixture-model mitochondrial filtering, reimplemented and switched off by default in steps 02 and 04.
- [CellBender](https://doi.org/10.1038/s41592-023-01943-7) (Fleming et al. 2023) and [SoupX](https://doi.org/10.1093/gigascience/giaa151) (Young and Behjati 2020). Ambient RNA, the reason the dual-TRB filter has a second-chain floor.

## Gene lists

- [Dissecting the multicellular ecosystem of metastatic melanoma](https://doi.org/10.1126/science.aad0501) (Tirosh et al. 2016). The S and G2M cell-cycle gene sets used in the lateral list.
- [Single-cell sequencing reveals dissociation-induced gene expression](https://doi.org/10.1038/nmeth.4437) (van den Brink et al. 2017). Why immediate-early and heat-shock genes are treated as lateral.

## Annotation

- [ComplexHeatmap reference book](https://jokergoo.github.io/ComplexHeatmap-reference/book/). The module and annotation heatmaps in steps 03, 05 and 05b are built on it; `column_km` is the block-count setting.
- [Azimuth references](https://azimuth.hubmapconsortium.org/). Marker tables per cell type and subtype. Useful for checking which markers cannot separate two states, and therefore belong in `panel` (display) rather than `sig` (scoring).

## Composition

- [Compositional analysis chapter](https://www.sc-best-practices.org/conditions/compositional.html), sc-best-practices. The clearest explanation of why proportions cannot be tested directly.
- [Microbiome datasets are compositional: and this is not optional](https://doi.org/10.3389/fmicb.2017.02224) (Gloor et al. 2017). Written for microbiome data, but the clearest short account of the centred log-ratio and how to read it.
- [propeller](https://doi.org/10.1093/bioinformatics/btac582) (Phipson et al. 2022) and [scCODA](https://doi.org/10.1038/s41467-021-27150-6) (Büttner et al. 2021). Two dedicated alternatives to the CLR and mixed model used in step 07.

## Differential expression

- [Confronting false discoveries in single-cell differential expression](https://doi.org/10.1038/s41467-021-25960-2) (Squair et al. 2021). The paper to read if pseudobulk still looks optional.
- [A practical solution to pseudoreplication bias in single-cell studies](https://doi.org/10.1038/s41467-021-21038-1) (Zimmerman et al. 2021) and [A balanced measure shows superior performance of pseudobulk methods](https://doi.org/10.1038/s41467-022-35519-4) (Murphy and Skene 2022). The debate over how to model donors, and why the donor is the unit.
- [DESeq2 paper](https://doi.org/10.1186/s13059-014-0550-8) (Love et al. 2014) and [vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). Design formulas, contrasts and shrinkage.
- [False discovery rates: a new deal](https://doi.org/10.1093/biostatistics/kxw041) (Stephens 2017). The `ashr` method used to shrink fold changes in step 08.
- [From reads to genes to pathways](https://doi.org/10.12688/f1000research.8987.2) (Chen, Lun and Smyth 2016). The edgeR workflow that introduced `filterByExpr`.
- [dream](https://doi.org/10.1093/bioinformatics/btaa687) (Hoffman and Roussos 2021). The mixed-model method used for the interaction family.

## Pathway enrichment

- [Gene set enrichment analysis](https://doi.org/10.1073/pnas.0506580102) (Subramanian et al. 2005). The original method, including the running enrichment score and the leading edge.
- [The Molecular Signatures Database Hallmark gene set collection](https://doi.org/10.1016/j.cels.2015.12.004) (Liberzon et al. 2015). How the 50 Hallmark sets were built to reduce redundancy.
- [fgsea](https://doi.org/10.1101/060012) (Korotkevich et al.) and its [vignette](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html). The fast algorithm under `clusterProfiler::GSEA`.
- [clusterProfiler 4.0](https://doi.org/10.1016/j.xinn.2021.100141) (Wu et al. 2021) and the [clusterProfiler book](https://yulab-smu.top/biomedical-knowledge-mining-book/).
- [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/). What is in Hallmark, GO BP, KEGG and Reactome.

## Clonal repertoire

- [scRepertoire documentation](https://www.borch.dev/uploads/screpertoire/). Vignettes for every function in steps 10 and 11.
- [rnabio.org clonality module](https://rnabio.org/module-08-scrna/0008/06/01/Clonality/). The workflow step 10 is structured after.
- [Mhanna et al. 2024](https://doi.org/10.1016/j.crmeth.2024.100738). Occurrence-based clone size bins, Renyi profiles and rarefaction curves for repertoire data.
- [Lineage tracking reveals dynamic relationships of T cells in colorectal cancer](https://doi.org/10.1038/s41586-018-0694-x) (Zhang et al. 2018). The origin of the Startrac expansion, migration and transition indices.
- [TCR beta allelic exclusion](https://doi.org/10.4049/jimmunol.1001158) (Brady et al. 2010) and [dual TCR alpha expression](https://doi.org/10.1126/science.8493531) (Padovan et al. 1993). Why beta, not alpha, is used to call doublets.
- [10x V(D)J algorithm documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-5p-vdj-algorithm). What the `filtered_contig_annotations.csv` columns mean.

## Tools

- [Seurat v5](https://satijalab.org/seurat/)
- [scanpy](https://scanpy.readthedocs.io/) and [anndata](https://anndata.readthedocs.io/)
- [reticulate](https://rstudio.github.io/reticulate/), how the R steps read `.h5ad` files
- [Apptainer](https://apptainer.org/docs/user/latest/) for the containers in [14](14_reproducibility.md)
- [Posit Package Manager](https://packagemanager.posit.co/), the dated CRAN snapshots used to pin the R image
