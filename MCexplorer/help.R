# =============================================================================
#  MCexplorer  |  help.R
#  The "How to use" manual of each page, shown in a window:
#    what     one line on the page
#    steps    3 to 5 cards: c(icon, title, text)
#    visual   optional drawing (HELP_VISUAL in global.R)
#    example  list(flow = c("Label|value", ...), result = "...")
#    manual   sections: list(icon, title, items = list(c(control, explanation), ...))
#             'control' is drawn like the control on the page; an empty control gives a plain paragraph
#  A little HTML (<b>, <i>) is allowed.
# =============================================================================

HELP <- list(
  home = list(
    title = "Home: opening a dataset",
    what = "Point the app at the output folder of the gene-module step; cell information and a UMAP are optional extras.",
    steps = list(
      c("folder-open", "Paste the folder", "Its full path. The files are found as soon as the path is typed."),
      c("file-circle-check", "Check the files", "Four module files are needed; several runs in one folder can be chosen from a list."),
      c("users", "Add cell information", "The annotated Seurat object or a table of cells: labels, samples, conditions."),
      c("arrow-right", "Open", "Large datasets take about a minute the first time, then seconds.")),
    example = list(flow = c("Folder|/sc/arion/.../MC_modules", "Found|4 module files", "Cell information|preCART_annotated.rds", "UMAP|compute"),
                   result = "Open dataset, then the Annotate page."),
    manual = list(
      list("folder-tree", "The files", list(
        c("genelist .rd", "The genes of every module (output of the gene-module step)."),
        c("metadata .rd", "Module scores of every metacell (the <i>_scaled</i> columns of metadf_z)."),
        c("counts .rd", "Counts per gene and metacell (dgcmtxcounts); used for genes, markers and mixed blocks."),
        c("HMlist .rd", "The pipeline heatmap; its column splits are the <b>pipeline blocks</b>."))),
      list("users", "Cell information", list(
        c("Seurat .rds", "Read in the background and summarised per metacell, so the app never keeps all cells in memory."),
        c("Table of cells", "Any .csv or .tsv with a metacell column and cell labels (broad_label, sample, condition...)."),
        c("", "Without it the app still works, but the Conditions page and the label bars are empty."))),
      list("gauge-high", "Speed", list(
        c("", "Slow steps are saved next to the data files (a small .mcexplorer file), so a second opening takes seconds."),
        c("", "Uploading through the browser works but is slow for large Seurat objects; a folder path is faster."),
        c("", "A folder can also be given when starting: MCX_DATA=/path Rscript -e 'shiny::runApp(\"MCexplorer\")'."))))),

  annotate = list(
    title = "Annotate",
    what = "Group metacells into blocks, look at what each block contains, and give it a name. Metacells are never split or changed: names are given to whole metacells.",
    steps = list(
      c("layer-group", "Choose blocks", "Pipeline blocks, new blocks (any number), or an existing label."),
      c("hand-pointer", "Hover and click", "Hovering shows a block's size, main label and name; a click selects it."),
      c("magnifying-glass", "Read the panel", "Labels with their share of cells, names already given, highest modules, top genes."),
      c("tag", "Name it", "Type and press Enter. Mixed? Split the block first.")),
    visual = "annotate",
    example = list(flow = c("Block|B22", "broad_label|84% Mono/Mac", "Modules|M12 C1QA C1QB", "Name|Resident macrophages"),
                   result = "B23 is only 55% Mono/Mac and flagged <i>mixed</i>: Check if mixed, Split into 2, then name each part."),
    manual = list(
      list("layer-group", "Blocks bar (above the heatmap)", list(
        c("Pipeline blocks", "The column splits of the pipeline heatmap (ComplexHeatmap column_km)."),
        c("New blocks", "All metacells re-grouped with the same method, into any number of blocks. The pipeline's number of blocks with seed 42 gives the pipeline blocks back exactly; more blocks give finer groups."),
        c("seed", "k-means starts from random centres; the seed fixes the result so it can be repeated."),
        c("An existing label", "Metacells grouped by e.g. their main broad_label, to check or refine an existing annotation."),
        c("Every metacell", "No blocks: each metacell is selected on its own."))),
      list("table-cells", "The heatmap", list(
        c("hover", "A box shows the block under the mouse with its size, main label and the name it has."),
        c("click", "Selects the block. Clicking it again deselects it."),
        c("drag", "Selects every block the drag touches."),
        c("Shift / Ctrl", "Held while clicking or dragging, adds to the selection instead of replacing it."),
        c("Esc / click outside", "Clears the selection (clicks on buttons and panels do not)."),
        c("double-click", "Zooms into the block under the mouse; double-clicking again zooms back out."),
        c("Zoom to selection", "Shows only the selected blocks across the full width, to look at them closely."))),
      list("sliders", "Left panel", list(
        c("Columns", "<b>Every metacell</b> draws one column per metacell; <b>Block averages</b> one column per block, easier with thousands of metacells."),
        c("Label bars", "Coloured bars above the heatmap: labels from the Seurat object, pipeline blocks, names given."),
        c("Modules shown", "All modules, or chosen ones. Their genes are listed under the heatmap, with a search box."),
        c("Size", "Height of the heatmap and the colour range of the module scores."))),
      list("tag", "Selection panel (right)", list(
        c("label shares", "Share of the selection's cells in each value of the main label, e.g. 84% Mono/Mac. Below it, the metacells counted by their own main label: a block of mixed metacells can look uniform on the label bar while its cells are split."),
        c("Names already given", "Shown when part of the selection is named already, so a name is not overwritten by accident."),
        c("Apply name", "Gives the typed name to every selected metacell. Names used before are suggested while typing."),
        c("Split", "Splits one block into parts by re-clustering only its own metacells; the parts are then named separately."),
        c("Merge", "Shift or Ctrl click two or more blocks that hold the same population, then Merge: they become one block (e.g. B3+B7). Only the grouping changes, never the metacells."),
        c("Undo last / Undo all", "In the Blocks bar: take back the last split or merge, or all of them."),
        c("Check if mixed", "Opens Mixed groups on Modules & markers with this block chosen."))),
      list("list-check", "Names given (right, below)", list(
        c("name", "Click to select its metacells; the pencil renames, the cross sets the name back to undefined."),
        c("Start from a label", "Copies e.g. broad_label as first names to correct block by block."),
        c("Load names", "Reads names saved earlier: the Excel file from Export, or a CSV with metacell and annotation columns."),
        c("", "Every change is saved automatically; Undo reverts the last one."))))),

  heat2 = list(
    title = "Heatmap & genes",
    what = "Expression of chosen genes in the same columns as the Annotate heatmap, to check that a block expresses the genes expected for its name.",
    steps = list(
      c("list", "Choose genes", "Top genes per block, a module, a signature, a gene set or a typed list."),
      c("bars", "Label bars", "Blocks, names and labels shown as coloured bars on top."),
      c("table-cells", "Module scores", "Optional: module scores under the genes."),
      c("palette", "Adjust", "Colours top right; scaling and layout on the left.")),
    example = list(flow = c("Genes|CD68 C1QA CD14 S100A8", "Bars|Annotate blocks, broad_label", "Module scores|M12, M7"),
                   result = "Blocks named Resident macrophages show C1QA and CD68 but not S100A8."),
    manual = list(
      list("list", "Genes", list(
        c("Top genes of each block", "The genes most specific to each block against all other metacells; a quick overview."),
        c("Most variable genes", "Genes that differ most between blocks."),
        c("Genes of a module", "Every gene of one module."),
        c("A signature or gene set", "A typed list (+) or a published gene set (book: Hallmark, GO, Reactome, KEGG, or a .gmt file)."),
        c("A typed gene list", "Any genes, separated by commas, spaces or new lines."))),
      list("table-cells", "Display", list(
        c("Columns", "Follow the blocks of the Annotate page; the link above the plot goes there."),
        c("Module scores", "Off by default. All or chosen modules drawn under the genes for comparison."),
        c("Scale each gene", "Every gene coloured by its own range, so weakly expressed genes show a pattern."),
        c("Column order from a file", "Draws metacells in a fixed order from an Excel file, e.g. to match a figure from a script."))))),

  plot = list(
    title = "Expression by group",
    what = "Modules, genes or signatures across groups of metacells, as box, violin, bar, strip, heatmap or dot plots, with optional statistics.",
    steps = list(
      c("layer-group", "Pick groups", "A label, the names given, or the Annotate blocks. Drag to reorder."),
      c("dna", "Pick what to plot", "One or more modules, genes or signatures."),
      c("chart-simple", "Pick a plot type", "Above the plot; hover a type for what it shows."),
      c("square-root-variable", "Statistics", "Switch on above the plot; settings open on the left.")),
    example = list(flow = c("Groups|broad_label", "Genes|IFNG, CXCL9", "Type|Violin", "Statistics|vs T cells"),
                   result = "Two panels with Wilcoxon p-values, Holm-adjusted."),
    manual = list(
      list("chart-simple", "Plot types", list(
        c("Box", "Median, quartiles and range."), c("Violin", "Shape of the distribution."), c("Bar", "Mean with standard error."),
        c("Strip", "Every metacell with a median line."), c("Heatmap", "Mean per group; one row per item."),
        c("Dot plot", "Colour: mean; size: share of metacells above zero."),
        c("", "Metrics such as the number of cells per metacell are single values, so only box, violin, bar and strip plots apply."))),
      list("dna", "Signatures", list(
        c("Mean z-score", "Each gene scaled across metacells, then averaged: simple, every gene counts equally."),
        c("Matched controls", "Average minus that of control genes with similar overall expression (Tirosh 2016, as Seurat AddModuleScore); removes differences in overall detection."),
        c("", "The method of the chosen signature is shown, and can be changed, right under the signature box."))),
      list("square-root-variable", "Statistics", list(
        c("Comparisons", "All pairs, each group against a reference, or chosen pairs."),
        c("Test", "Wilcoxon (rank-based, the usual choice) or Welch's t-test."),
        c("Correction", "Holm (default), Benjamini-Hochberg (FDR), Bonferroni or none."),
        c("Heatmap / dot plot", "In every row, each column is tested against the reference column; the result is written on the tile."),
        c("", "Each point is one metacell, so these p-values rank differences rather than prove them; condition effects belong on the Conditions page."))))),

  cond = list(
    title = "Conditions",
    what = "Composition or expression per sample, compared between conditions with samples as replicates.",
    steps = list(
      c("chart-pie", "Composition or expression", "The two large buttons at the top."),
      c("vials", "Conditions", "Condition column, which conditions, and one point per sample."),
      c("layer-group", "Cell groups", "Which groups to plot, or which cells to average."),
      c("ranking-star", "Ranked table", "Every group tested; biggest differences first.")),
    example = list(flow = c("Condition|group: EC, CTRL, UC, NV", "One point per|sample", "Groups|Names given"),
                   result = "Inflammatory monocytes rank first; plotting them shows a higher share in EC."),
    manual = list(
      list("chart-pie", "Two questions", list(
        c("Cell composition", "Share of each sample's cells that fall in a group. All cells of a sample are the denominator, so shares of different groups compare directly."),
        c("Expression", "A module, gene or signature averaged over the chosen cells of each sample, each metacell weighted by its number of cells."))),
      list("sliders", "Left panel", list(
        c("Condition column", "The column that holds each sample's condition."),
        c("One point per", "The replicate: usually the sample or patient."),
        c("Groups from", "Names given on Annotate, a label, or the blocks."),
        c("Sample filter", "Keeps only samples with chosen values of another column, e.g. only colon biopsies."))),
      list("table", "Tables (under the plot)", list(
        c("Ranked", "Every group (or module) tested across the conditions shown, biggest differences first. Tick rows and plot them."),
        c("Values per sample", "The numbers behind the plot."),
        c("Statistics", "Every comparison drawn on the plot."))))),

  genes = list(
    title = "Modules & markers",
    what = "Four gene lookups. The left panel changes with the button chosen.",
    steps = list(
      c("layer-group", "Modules", "Genes of a module and where it is high."),
      c("shuffle", "Mixed groups", "Genes varying inside a block: one cell type or two?"),
      c("star", "Marker genes", "Genes that set a group apart."),
      c("magnifying-glass", "Gene lookup", "Module of a gene and where it is expressed.")),
    example = list(flow = c("Mixed groups|B23", "Top genes|S100A8, C1QA, FCN1"),
                   result = "Monocyte and macrophage genes vary inside B23: split it on Annotate."),
    manual = list(
      list("shuffle", "Mixed groups", list(
        c("what it tests", "For each gene, whether its counts differ between the metacells of a block more than their sizes explain (chi-square on counts)."),
        c("Chi-square", "Higher means more variation between the block's metacells."),
        c("reading it", "Many strong genes from two programmes (e.g. S100A8 and C1QA) mean the block mixes cell states: split it. Few or weak genes mean it is homogeneous."),
        c("", "Reached directly from Annotate with <i>Check if mixed</i> under a selected block."))),
      list("star", "Marker genes", list(
        c("Compared with", "All other metacells (what defines a group) or chosen groups (what separates two similar groups)."),
        c("log2 fold change", "Mean expression in the group against the comparison, the usual metacell convention."))),
      list("layer-group", "Modules and gene lookup", list(
        c("Module", "Type a number or any gene in it; the plot shows where the module is high and its genes are listed under the plot."),
        c("Gene", "Any measured gene, also genes in no module."))))),

  umap = list(
    title = "UMAP",
    what = "Metacells in two dimensions, coloured by anything known about them.",
    steps = list(
      c("palette", "Colour by", "A label, the names given, a module, a gene or a signature."),
      c("font", "Names", "Group names at the centre of each group."),
      c("ruler-combined", "Layout", "Recompute with other settings.")),
    manual = list(
      list("braille", "Layout", list(
        c("Neighbours", "More neighbours keep the large-scale structure; fewer keep local detail."),
        c("Min. distance", "How tightly points are packed."),
        c("", "A UMAP .csv chosen on Home is used as is; otherwise one is computed from the module scores with uwot."))))),

  export = list(
    title = "Export",
    what = "Take out the names, the figures and tables, or the whole session.",
    steps = list(
      c("file-excel", "Annotation", "Excel with a summary sheet and a Metacells sheet for the Seurat object."),
      c("images", "Figures and tables", "Tick what goes in one zip."),
      c("floppy-disk", "Session", "Save and resume later.")),
    example = list(flow = c("Excel|Metacells sheet", "R|match() on the metacell column"),
                   result = "obj$new_annotation &lt;- m$annotation[match(obj$metacell, m[[1]])]"),
    manual = list(
      list("file-excel", "The annotation file", list(
        c("Annotation sheet", "One row per name: metacells, cells, main labels, top modules, top genes, share of each condition's cells."),
        c("Metacells sheet", "One row per metacell with its name; the first column is named like the metacell column of the Seurat object."))),
      list("images", "Figures and tables", list(
        c("Figures", "Every figure opened in this session, at the size set in the Figure bar under it, as PDF, PNG (300 dpi) or SVG."),
        c("Tables", "Saved as sheets of one Excel file."))))))
