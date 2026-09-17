# Gene lists

MC2 takes three separate gene lists and they do three different things. Mixing
them up is the most common way to end up with metacells that look fine and are
not.

## Excluded

Genes dropped from the object entirely. Not counted, not used for similarity,
not present downstream.

For things that are not the biology being measured:

- **Mitochondrial genes** (`^MT-`). High mito fraction is a QC signal, not a
  cell state, and leaving these in lets dying cells cluster together across
  every lineage.
- **`MALAT1`, `NEAT1`.** Nuclear lncRNAs whose count is dominated by how much
  nucleus made it into the droplet.
- **Haemoglobin** (`HBB`, `HBA1`, `HBA2`, `HBD`, `HBM`). Ambient from lysed red
  cells, present at some level in nearly every PBMC droplet.
- **`MTRNR2L*`.** Nuclear-encoded humanin-like pseudogenes. Sequence-similar to
  the mitochondrial `MT-RNR2` and heavily expressed in stressed cells.

> [!WARNING]
> Do not use the pattern `MT1.*` for mitochondrial genes. It matches the
> metallothioneins `MT1A`, `MT1X`, `MT2A`, which are real stress-response
> genes. Use `^MT-` with the hyphen.

Excluding a gene also changes `total_umis`, because the total is computed after
exclusion. That number is what `target_metacell_umis` is compared against, so
the post-exclusion total is the one that matters when setting it.

## Lateral

Genes still counted in totals and still used for outlier detection, but not
used to compute cell-cell similarity.

This is the subtle one. A lateral gene is real biology that should not define
cell identity. Cell cycle is the archetype: a proliferating CD8 T cell and a
proliferating monocyte share a strong cell-cycle programme, and if those genes
drive the graph the result is a cycling metacell containing both. Marking them
lateral keeps the cycling signal measurable in the output while stopping it
pulling two lineages together.

The list here:

| Group | Why lateral |
|:--|:--|
| Cell cycle S and G2M ([Tirosh et al. 2016](https://doi.org/10.1126/science.aad0501)) | Shared across every proliferating lineage |
| Stress and immediate-early (`FOS`, `JUN`, `HSPA*`, `DNAJB1`) | Dissociation artefact ([van den Brink et al. 2017](https://doi.org/10.1038/nmeth.4437), [O'Flanagan et al. 2019](https://doi.org/10.1186/s13059-019-1830-0)) |
| Sex-linked (`XIST`, `RPS4Y1`, `DDX3Y`) | Splits every cell type by donor sex |
| Platelet and megakaryocyte (`PPBP`, `PF4`) | Ambient from platelet contamination |
| HLA class I and II | Varies with activation and donor genotype, not identity |
| Ribosomal (`^RPS`, `^RPL`) | Tracks library size and translational state |
| TCR and BCR V/J segments | See below |
| Immunoglobulin V/J | Same argument as TCR |

### The V/J segments

`TRAV`, `TRBV`, `TRAJ`, `TRBJ` and their BCR equivalents are the variable
segments of the receptor. Every cell in a clone carries the same ones. Left in
the similarity computation, the graph starts grouping cells by clone rather
than by state, and a metacell comes out as one expanded clone wearing a cell
state's clothes.

This matters more on a T cell subset than it did globally, not less. Globally,
between-lineage variance swamps it. Inside a single lineage, clonal variance is
one of the largest signals present, and it is exactly the wrong one.

The constant regions `TRAC`, `TRBC1`, `TRBC2` are kept. They are shared across
clones and are useful lineage markers. The same logic drives
[scRepertoire's `quietTCRgenes`](https://doi.org/10.1093/bfgp/elac049).

## Noisy

Genes given extra slack in outlier detection because they are bursty.

Orthogonal to lateral. A noisy gene can still drive similarity; it is just not
allowed to make a cell look like an outlier on its own. Immunoglobulin constant
regions are the clearest case: a plasma cell's `IGHG1` count can sit two orders
of magnitude above its neighbours' without that cell being anything other than
a plasma cell.

Here the noisy list is the stress genes plus the immunoglobulin constant
regions plus a few high-dynamic-range HLA genes. The overlap with the lateral
list is intentional; the two flags do not conflict.

## The check that catches mistakes

After feature selection, steps 02 and 04 both run:

```python
leaked = clean.var_names[clean.var["selected_gene"] & clean.var["lateral_gene"]]
```

This should be empty. A lateral gene in the selection means the lateral marking
did not take effect, usually because it was applied to the wrong object or
after selection had already run. A leaked cell-cycle gene visibly degrades the
metacells. Add it, re-mark, rerun.

Step 04 adds a stricter version for the subset case: selected genes matching
`TR[ABGD][VJD]`. Anything printed there on a T cell subset is a stop, not a
warning.

> [!TIP]
> Marking a gene lateral is cheap and reversible. Excluding it is not, because
> it changes every total. When unsure which list something belongs in, lateral
> is almost always the safer answer.
