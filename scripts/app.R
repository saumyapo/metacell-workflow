# =============================================================================
#  MCexplorer: explore metacell gene-module results
# =============================================================================
#
#  HOW TO START
#    RStudio:   paste the folder that holds your metacell files into DATA_PATH below, then click "Run App".
#    Terminal:  MCX_DATA=/full/path/to/folder Rscript -e 'shiny::runApp("app.R", port = 8888, launch.browser = FALSE)'
#    Or leave DATA_PATH empty and open / upload the files from the app's Data page.
#
#  WHAT THE FOLDER (OR UPLOAD) NEEDS
#    Required: the 4 files written by the metacell gene-module step
#      ...genelist...rd   gene modules              (object: genelist)
#      ...metadata...rd   module scores per metacell (object: metadf_z)
#      ...counts...rd     gene expression            (object: dgcmtxcounts)
#      ...HMlist...rd     the module heatmap         (object: ht)
#    Optional:
#      cell information: a Seurat object (.RDS) or a table of cells (.csv / .tsv) with a column of metacell names.
#                        Gives cell types, samples, conditions and cell counts. Read once in the background, then cached.
#      metacell UMAP:    a .csv with metacell names and two coordinate columns. If missing, the app computes one.
#
#  The first opening of a dataset is the slow one; results are cached next to the data (or in a temp folder).
# =============================================================================

DATA_PATH <- ""

if (!nzchar(DATA_PATH)) DATA_PATH <- Sys.getenv("MCX_DATA", "")

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(Matrix)
  library(DT)
})

options(shiny.maxRequestSize = 20 * 1024^3)

APP_VERSION  <- "5.0"
EPS          <- 1e-5
CACHE_VER    <- 1
AUTOSAVE_DIR <- tools::R_user_dir("MCexplorer", "data")

# gene list used by the original app (v3.8), offered as an option only
ORIGINAL_LIST <- c("CD2","CD7","IL7R","CD19","CD14","NCAM1","CD3G","CD3D","CD3E","CD4","CD8A","CD8B","TRAC","TRBC1","TRBC2",
  "TRGC1","TRGC2","TRDC","IL2RA","FOXP3","TIGIT","CTLA4","IL10","IL18R1","IL18RAP","SLC4A10","KLRB1","ZBTB16","IL23R",
  "TRAV1-2","ZNF683","CCR7","SELL","TCF7","LEF1","KLF2","S1PR1","DUSP2","CCL5","CD27","CD28","FAS","PDCD1","B3GAT1",
  "HLA-DRA","CD38","FOS","JUN","NR4A1","NR4A2","NR4A3","CD69","IL2","TBX21","IFNG","TNF","CSF2","XCL1","XCL2","GATA3",
  "IL4","IL5","IL13","RORC","IL17A","IL21","IL22","FCGR3A","GNLY","PRF1","NCR1","EOMES","NKG7","KLRC1","KLRK1",
  "KIR2DL1","GZMA","GZMB","GZMH","CD200","BTLA","CCR4")

PAL <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#86BCB6",
         "#F1CE63","#D37295","#A0CBE8","#FFBE7D","#8CD17D","#B6992D","#499894","#FABFD2","#D4A6C8","#2F4B7C",
         "#A05195","#665191","#D45087","#79706E")
GREY <- "#C9CDD4"
NO_LABEL <- c("undefined", "Undefined", "Unassigned", "Outlier", "NA")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ============================================================== small helpers

mcx_log <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
fmt_secs <- function(s) { s <- as.numeric(s); if (s < 60) sprintf("%.0f s", s) else sprintf("%d min %02d s", as.integer(s %/% 60), as.integer(round(s %% 60))) }
fmt_size <- function(b) format(structure(b, class = "object_size"), units = "auto")
fmt_n    <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
pretty_name <- function(x) gsub("_", " ", x)

parse_list <- function(txt) {
  if (is.null(txt) || !nzchar(trimws(txt))) return(character(0))
  s <- unlist(strsplit(gsub("^\\s*c\\(|\\)\\s*$", "", trimws(txt)), "[,;\n\t ]+"))
  s <- trimws(gsub("[\"']", "", s))
  unique(s[nzchar(s)])
}

match_genes <- function(g, universe) {
  hit <- ifelse(g %in% universe, g, universe[match(toupper(g), toupper(universe))])
  list(found = unique(hit[!is.na(hit)]), missing = g[is.na(hit)])
}

mixed_sort <- function(x) {
  x <- unique(as.character(x[!is.na(x)]))
  n <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", x)))
  x[order(gsub("[0-9.]", "", x), n, x)]
}

file_key <- function(paths) paste(basename(paths), file.size(paths), as.numeric(file.mtime(paths)), collapse = "|")

chunks <- function(n, size = 1000) split(seq_len(n), ceiling(seq_len(n) / size))

# loads one .rd file and returns the object, by expected name or the only object present
load_rd <- function(path, expected) {
  e <- new.env()
  nm <- load(path, envir = e)
  if (expected %in% nm) return(get(expected, envir = e))
  if (length(nm) == 1) return(get(nm, envir = e))
  stop(basename(path), " does not contain '", expected, "'")
}

# role of an object loaded from an .rd file
rd_role <- function(obj) {
  if (inherits(obj, c("HeatmapList", "Heatmap"))) return("heatmap")
  if (inherits(obj, "Matrix") || (is.matrix(obj) && is.numeric(obj))) return("counts")
  if (is.data.frame(obj) && any(grepl("_scaled$", colnames(obj)))) return("scores")
  if (is.data.frame(obj) && ncol(obj) >= 2 && is.character(obj[[2]])) return("genelist")
  NA_character_
}

ROLE_KEYS  <- c(genelist = "genelist", scores = "metadata", counts = "counts", heatmap = "hmlist")
ROLE_NAMES <- c(genelist = "Gene modules", scores = "Module scores", counts = "Gene expression", heatmap = "Module heatmap")
ROLE_OBJ   <- c(genelist = "genelist", scores = "metadf_z", counts = "dgcmtxcounts", heatmap = "ht")

# finds complete sets of the 4 files in a folder, grouped by their shared name
scan_folder <- function(dir) {
  f  <- list.files(dir, full.names = TRUE)
  rd <- f[grepl("\\.(rd|rda|rdata)$", f, ignore.case = TRUE)]
  role <- vapply(basename(rd), function(x) { k <- names(ROLE_KEYS)[vapply(ROLE_KEYS, grepl, TRUE, x = tolower(x))]; if (length(k) == 1) k else NA_character_ }, "")
  stem <- tolower(sub("\\.(rd|rda|rdata)$", "", basename(rd), ignore.case = TRUE))
  for (k in ROLE_KEYS) stem <- sub(k, "", stem, fixed = TRUE)
  stem <- gsub("_+", "_", gsub("^_|_$", "", stem))
  ok <- !is.na(role)
  sets <- lapply(split(seq_along(rd)[ok], stem[ok]), function(i) setNames(rd[i], role[i]))
  complete <- Filter(function(s) all(names(ROLE_KEYS) %in% names(s)), sets)
  incomplete <- Filter(function(s) !all(names(ROLE_KEYS) %in% names(s)), sets)
  others <- f[!grepl("\\.(rd|rda|rdata)$", f, ignore.case = TRUE) & !grepl("^\\.mcexplorer", basename(f)) & !dir.exists(f)]
  cell_cand <- others[grepl("\\.(rds|csv|tsv|txt)(\\.gz)?$", others, ignore.case = TRUE) & !grepl("umap", basename(others), ignore.case = TRUE)]
  umap_cand <- others[grepl("umap", basename(others), ignore.case = TRUE) & grepl("\\.(csv|tsv|txt)$", others, ignore.case = TRUE)]
  list(dir = dir, sets = complete, incomplete = incomplete, cell_cand = cell_cand, umap_cand = umap_cand)
}

guess_cell_file <- function(cand, stem) {
  if (!length(cand)) return("")
  a <- cand[grepl("annot", basename(cand), ignore.case = TRUE)]
  if (!length(a)) return("")
  tok <- strsplit(stem, "_")[[1]]
  a[which.max(vapply(a, function(x) sum(vapply(tok, grepl, TRUE, x = tolower(basename(x)), fixed = TRUE)), 0))]
}

# ============================================================== data preparation

reporter <- function(n, title) {
  i <- 0; t0 <- Sys.time(); tl <- Sys.time()
  mcx_log("Preparing ", title)
  function(label) {
    if (i > 0) mcx_log(sprintf("      done (%s)", fmt_secs(difftime(Sys.time(), tl, units = "secs"))))
    i <<- i + 1; tl <<- Sys.time()
    mcx_log(sprintf("[%d/%d] %s", i, n, label))
    if (!is.null(getDefaultReactiveDomain())) incProgress(1 / (n + 1), detail = label)
    invisible(t0)
  }
}

# per-block sums used for markers, default gene lists and variability, one pass over the matrix
block_stats <- function(counts, libsize, block, nb) {
  n <- ncol(counts); g <- nrow(counts)
  B <- sparseMatrix(i = seq_len(n), j = block, x = 1, dims = c(n, nb))
  S <- Fr <- matrix(0, g, nb); s1 <- s2 <- numeric(g)
  for (idx in chunks(n)) {
    m <- as(counts[, idx, drop = FALSE] %*% Diagonal(x = 1 / libsize[idx]), "CsparseMatrix")
    Bi <- B[idx, , drop = FALSE]
    Fr <- Fr + as.matrix(m %*% Bi)
    m@x <- log2(m@x + EPS) - log2(EPS)
    S <- S + as.matrix(m %*% Bi)
    s1 <- s1 + Matrix::rowSums(m); m@x <- m@x^2; s2 <- s2 + Matrix::rowSums(m)
  }
  sizes <- tabulate(block, nb)
  dimnames(S) <- dimnames(Fr) <- list(rownames(counts), seq_len(nb))
  list(S = S, frac = sweep(Fr, 2, sizes, "/"), sizes = sizes,
       gene_frac = rowSums(Fr) / n, gene_var = setNames((s2 - s1^2 / n) / max(n - 1, 1), rownames(counts)))
}

# genes ranked by log2 fold change of mean log2(fraction + 1e-5) between blocks A and B (B = all others if NULL)
rank_genes <- function(bs, A, B = NULL, min_frac = EPS) {
  sz <- bs$sizes; n <- sum(sz)
  SA <- rowSums(bs$S[, A, drop = FALSE]); nA <- sum(sz[A])
  if (is.null(B) || !length(B)) { SB <- rowSums(bs$S) - SA; nB <- n - nA }
  else { SB <- rowSums(bs$S[, B, drop = FALSE]); nB <- sum(sz[B]) }
  fa <- as.vector(bs$frac[, A, drop = FALSE] %*% sz[A]) / nA
  out <- data.frame(gene = rownames(bs$S), log2FC = SA / nA - SB / max(nB, 1), mean_frac = fa, stringsAsFactors = FALSE)
  out <- out[out$mean_frac >= min_frac, ]
  out[order(-out$log2FC), ]
}

prepare_dataset <- function(src, name, cache_file = NULL) {
  step <- reporter(6, name)
  t0 <- step("Reading gene modules")
  genelist <- if (is.character(src$genelist)) load_rd(src$genelist, "genelist") else src$genelist
  step("Reading module scores")
  metadf_z <- if (is.character(src$scores)) load_rd(src$scores, "metadf_z") else src$scores
  step("Reading the saved heatmap")
  ht <- if (is.character(src$heatmap)) load_rd(src$heatmap, "ht") else src$heatmap
  step(paste0("Reading gene expression", if (is.character(src$counts)) paste0(" (", fmt_size(file.size(src$counts)), ", the largest file)") else ""))
  counts <- if (is.character(src$counts)) load_rd(src$counts, "dgcmtxcounts") else src$counts
  counts <- as(counts, "CsparseMatrix")

  step("Checking that the files fit together")
  genelist <- as.data.frame(genelist); colnames(genelist)[1:2] <- c("Mods", "Gens")
  zcols <- grep("_scaled$", colnames(metadf_z), value = TRUE)
  if (!length(zcols)) stop("The module scores file has no '<module>_scaled' columns.")
  Z <- as.matrix(metadf_z[, zcols, drop = FALSE])
  if (!identical(rownames(Z), colnames(counts))) {
    common <- intersect(rownames(Z), colnames(counts))
    if (length(common) < nrow(Z)) stop("Metacell names differ between the module scores and gene expression files.")
    counts <- counts[, rownames(Z)]
  }
  if (inherits(ht, "Heatmap")) ht <- draw(ht)
  row_ord <- as.integer(unlist(row_order(ht)))
  col_ord <- lapply(column_order(ht), as.integer); names(col_ord) <- NULL
  if (length(unlist(col_ord)) != nrow(Z) || length(row_ord) != ncol(Z)) stop("The saved heatmap does not match the module scores (different number of metacells or modules).")

  mcs <- rownames(Z); nb <- length(col_ord)
  block <- integer(length(mcs)); for (k in seq_len(nb)) block[col_ord[[k]]] <- k; names(block) <- mcs
  num <- suppressWarnings(as.integer(sub("^.*_(\\d+)_scaled$", "\\1", zcols))); if (anyNA(num)) num <- seq_along(zcols)
  key <- sub("_scaled$", "", zcols)
  gtxt <- genelist$Gens[match(key, genelist$Mods)]; gtxt[is.na(gtxt)] <- ""
  mods <- data.frame(id = zcols, num = num, name = paste("Module", num), genes = gtxt,
                     n_genes = lengths(strsplit(gtxt, ",")), stringsAsFactors = FALSE)
  libsize <- Matrix::colSums(counts)
  allsum <- setNames(unlist(lapply(chunks(ncol(counts)), function(i) { s <- counts[, i, drop = FALSE]; s@x <- log1p(s@x); Matrix::colSums(s) }), use.names = FALSE), mcs)

  ckey <- if (is.character(src$counts)) file_key(unlist(src[c("genelist", "scores", "counts", "heatmap")])) else NULL
  cache <- cache_read(cache_file, ckey)
  if (!is.null(cache$bs)) { step("Summarising blocks (from cache)"); bs <- cache$bs }
  else { step("Summarising blocks (once; cached for next time)"); bs <- block_stats(counts, libsize, block, nb); cache_write(cache_file, ckey, "bs", bs) }
  mcx_log(sprintf("      done. Ready in %s: %s metacells, %d modules, %d blocks", fmt_secs(difftime(Sys.time(), t0, units = "secs")), fmt_n(length(mcs)), ncol(Z), nb))

  block_z <- rowsum(Z, block) / tabulate(block, nb)
  list(name = name, src = lapply(src, function(x) if (is.character(x)) x else "(uploaded)"), opened = Sys.time(),
       genelist = genelist, Z = Z, counts = counts, ht = ht, row_ord = row_ord, col_ord = col_ord,
       mcs = mcs, nb = nb, block = block, mc_order_ht = mcs[unlist(col_ord)], mods = mods,
       gene2mod = setNames(rep(mods$name, mods$n_genes), unlist(strsplit(gtxt, ","))),
       libsize = libsize, allsum = allsum, nfeat = diff(counts@p), bs = bs, block_z = block_z,
       block_top = lapply(seq_len(nb), function(k) head(rank_genes(bs, k, min_frac = 5e-5)$gene, 20)),
       cache_file = cache_file, cache_key = ckey, cache = cache)
}

cache_read <- function(file, key) {
  if (is.null(file) || is.null(key) || !file.exists(file)) return(list())
  x <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(x) || !identical(x$key, key) || !identical(x$ver, CACHE_VER)) return(list())
  x
}
cache_write <- function(file, key, field, value) {
  if (is.null(file) || is.null(key)) return(invisible())
  x <- cache_read(file, key); x$key <- key; x$ver <- CACHE_VER; x[[field]] <- value
  tryCatch(saveRDS(x, file), error = function(e) mcx_log("Could not write cache: ", conditionMessage(e)))
}
cache_path <- function(dir, stem) {
  d <- if (file.access(dir, 2) == 0) dir else tempdir()
  paste0(d, "/.mcexplorer_", gsub("[^A-Za-z0-9_]+", "_", stem), ".rds")
}

# Self-contained (runs in a background R process): reads cell information and summarises it per metacell.
summarise_cells <- function(path, mcs, max_levels = 200) {
  low <- tolower(sub("\\.gz$", "", path, ignore.case = TRUE))
  if (grepl("\\.(csv|tsv|txt)$", low)) {
    md <- if (requireNamespace("data.table", quietly = TRUE)) as.data.frame(data.table::fread(path))
          else utils::read.table(path, header = TRUE, sep = if (grepl("\\.csv$", low)) "," else "\t", check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    obj <- readRDS(path)
    md <- if (isS4(obj) && methods::.hasSlot(obj, "meta.data")) obj@meta.data
          else if (is.data.frame(obj)) obj
          else if (is.list(obj)) { d <- Filter(is.data.frame, obj); if (!length(d)) stop("No table of cells found in this file."); d[[which.max(vapply(d, nrow, 1L))]] }
          else stop("This file is not a Seurat object or a table.")
    rm(obj); invisible(gc())
  }
  md <- as.data.frame(md, stringsAsFactors = FALSE)
  cand <- names(md)[vapply(md, function(v) is.character(v) || is.factor(v), TRUE)]
  hits <- vapply(cand, function(c) sum(unique(as.character(md[[c]])) %in% mcs), 0)
  rn_hits <- sum(unique(rownames(md)) %in% mcs)
  if (rn_hits > max(c(hits, 0))) { md$..mc <- rownames(md); mc_col <- "..mc" }
  else {
    if (!length(hits) || max(hits) < 0.5 * length(mcs)) stop("No column in this file holds the metacell names used by the heatmap.")
    mc_col <- cand[which.max(hits)]
  }
  mc <- as.character(md[[mc_col]]); keep <- !is.na(mc) & mc %in% mcs
  md <- md[keep, , drop = FALSE]; mcf <- factor(mc[keep], levels = mcs)
  n_cells <- tabulate(mcf, length(mcs)); cell_level <- max(n_cells) > 1
  labels <- list(); comp <- list(); comp_major <- list(); comp_purity <- list(); nums <- list(); lv <- list(); total_umis <- NULL
  for (c in setdiff(names(md), c(mc_col, "..mc"))) {
    v <- md[[c]]
    if (is.numeric(v) && !is.factor(v)) {
      if (!cell_level || length(unique(v[seq_len(min(5000, length(v)))])) <= 1) next
      s <- as.vector(rowsum(as.numeric(v), mcf, reorder = TRUE, na.rm = TRUE)[as.character(mcs), 1])
      s[is.na(s)] <- 0
      if (is.null(total_umis) && grepl("^(ncount_rna|total_umis?)$", c, ignore.case = TRUE)) total_umis <- s
      nums[[c]] <- s / pmax(n_cells, 1)
      next
    }
    if (is.logical(v)) v <- as.character(v)
    u <- if (is.factor(v)) levels(droplevels(v)) else unique(v[!is.na(v)])
    if (length(u) < 2 || length(u) > max_levels) next
    f <- if (is.factor(v)) droplevels(v) else factor(v, levels = sort(u))
    tab <- unclass(table(mcf, f)); tab <- matrix(tab, nrow = length(mcs), dimnames = list(mcs, levels(f)))
    rs <- rowSums(tab); mx <- do.call(pmax, as.data.frame(tab))
    major <- colnames(tab)[max.col(tab, "first")]; major[rs == 0] <- NA
    lv[[c]] <- levels(f)
    if (!cell_level || sum(mx) / max(sum(rs), 1) >= 0.95) labels[[c]] <- major
    else { comp[[c]] <- tab; comp_major[[c]] <- major; comp_purity[[c]] <- round(mx / pmax(rs, 1), 3) }
  }
  list(file = path, mc_col = mc_col, cell_level = cell_level, n_total = nrow(md),
       n_cells = if (cell_level) n_cells else NULL, total_umis = total_umis,
       labels = labels, comp = comp, comp_major = comp_major, comp_purity = comp_purity, nums = nums, levels = lv)
}

read_umap <- function(path, mcs) {
  u <- if (grepl("\\.csv$", path, ignore.case = TRUE)) utils::read.csv(path, check.names = FALSE) else utils::read.delim(path, check.names = FALSE)
  idc <- which(vapply(u, function(v) sum(as.character(v) %in% mcs), 0) > 0.5 * length(mcs))
  ids <- if (length(idc)) as.character(u[[idc[1]]]) else rownames(u)
  num <- which(vapply(u, is.numeric, TRUE)); num <- setdiff(num, idc)
  if (length(num) < 2 || sum(ids %in% mcs) < 0.5 * length(mcs)) stop("Could not find metacell names and two coordinate columns in this file.")
  xy <- data.frame(x = u[[num[1]]], y = u[[num[2]]], row.names = ids)
  xy[intersect(rownames(xy), mcs), , drop = FALSE]
}

# ============================================================== plotting helpers

theme_mcx <- function(base = 13) {
  theme_minimal(base_size = base) + theme(
    plot.title = element_text(face = "bold", size = rel(1.1), colour = "#1D2129", margin = margin(b = 4)),
    plot.subtitle = element_text(colour = "#5B6170", size = rel(0.82), margin = margin(b = 8)),
    plot.caption = element_text(colour = "#8A90A0", size = rel(0.72)),
    axis.title = element_text(colour = "#3A3F4B", size = rel(0.88)),
    axis.text = element_text(colour = "#3A3F4B", size = rel(0.85)),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#ECEEF2", linewidth = 0.4),
    axis.line.x = element_line(colour = "#9AA0AC", linewidth = 0.4),
    legend.title = element_text(size = rel(0.82), face = "bold", colour = "#3A3F4B"),
    legend.text = element_text(size = rel(0.8)), legend.key.size = unit(0.9, "lines"),
    plot.margin = margin(12, 16, 10, 10))
}

save_plot <- function(file, draw, w, h, fmt) {
  w <- max(as.numeric(w %||% 8), 1); h <- max(as.numeric(h %||% 6), 1)
  if (fmt == "png") grDevices::png(file, width = w, height = h, units = "in", res = 300, type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
  else if (fmt == "svg") { if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite(file, width = w, height = h) else grDevices::svg(file, width = w, height = h) }
  else grDevices::pdf(file, width = w, height = h, useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  draw()
}

# ============================================================== UI helpers

tip <- function(text) tooltip(span(class = "mcx-tip", icon("circle-question")), text, placement = "right")
lab <- function(text, help = NULL) if (is.null(help)) text else tagList(text, tip(help))
section <- function(title, ..., help = NULL) div(class = "mcx-section", div(class = "mcx-section-title", title, if (!is.null(help)) tip(help)), ...)
hint <- function(...) div(class = "mcx-hint", ...)
seg <- function(...) div(class = "mcx-seg", ...)
page_head <- function(title, sub) div(class = "mcx-head", h2(title), p(sub))
empty_state <- function(ic, title, text) div(class = "mcx-empty", div(class = "i", icon(ic)), h5(title), p(text))
needs_data <- function(...) tagList(
  conditionalPanel("!output.has_data", empty_state("folder-open", "No dataset open yet", "Open a folder or upload files on the Data page.")),
  conditionalPanel("output.has_data", ...))
tab_body <- function(title, sub, sidebar_items, ...) div(class = "mcx-page", page_head(title, sub),
  needs_data(layout_sidebar(sidebar = sidebar(width = 340, open = "desktop", sidebar_items), border = TRUE, fillable = FALSE, ...)))

dl_bar <- function(id, w = 8, h = 6) {
  ns <- NS(id)
  div(class = "dl-bar",
      span(class = "dl-title", icon("download"), "Save figure"),
      selectInput(ns("fmt"), NULL, c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg"), width = "84px"),
      numericInput(ns("w"), NULL, w, min = 1, max = 80, step = 0.5, width = "72px"), span("x"),
      numericInput(ns("h"), NULL, h, min = 1, max = 80, step = 0.5, width = "72px"), span("inches"),
      downloadButton(ns("go"), "Download", class = "btn-sm btn-primary"))
}
dl_server <- function(id, draw, fname) moduleServer(id, function(input, output, session) {
  output$go <- downloadHandler(filename = function() paste0(fname(), ".", input$fmt),
                               content = function(file) save_plot(file, draw(), input$w, input$h, input$fmt))
})
plot_card <- function(plot_id, dl_id, w = 8, h = 6, ...) card(full_screen = TRUE, class = "mcx-plot-card",
  card_body(..., plotOutput(plot_id, height = "auto")), card_footer(dl_bar(dl_id, w, h)))

# ============================================================== CSS / JS

CSS <- "
:root { --mcx-bg:#F4F5F7; --mcx-ink:#1D2129; --mcx-muted:#6B7280; --mcx-line:#E5E7EB; --mcx-accent:#3557D6; --mcx-accent-soft:#EEF2FF; }
body { background:var(--mcx-bg); color:var(--mcx-ink); }
.navbar { background:#fff !important; border-bottom:1px solid var(--mcx-line); box-shadow:0 1px 2px rgba(16,24,40,.04); padding-top:0; padding-bottom:0; }
.navbar-brand { font-weight:700; letter-spacing:-.01em; color:var(--mcx-ink) !important; display:flex; align-items:center; gap:8px; }
.mcx-logo { width:22px; height:22px; border-radius:6px; background:linear-gradient(135deg,#3557D6 0%,#7C3AED 55%,#E15759 100%); display:inline-block; }
.navbar .nav-link { color:#4B5563 !important; font-weight:500; font-size:.88rem; padding:1rem .48rem !important; border-bottom:2px solid transparent; white-space:nowrap; }
.navbar > .container-fluid { padding-left:18px; padding-right:14px; }
.navbar-collapse { flex-wrap:nowrap; } .navbar-nav { flex-wrap:nowrap; }
.mcx-navstatus { display:flex; align-items:center; gap:4px; } .mcx-navstatus .mcx-chip { max-width:190px; overflow:hidden; text-overflow:ellipsis; margin-left:2px; }
.navbar .nav-link:hover { color:var(--mcx-ink) !important; }
.navbar .nav-link.active { color:var(--mcx-accent) !important; border-bottom-color:var(--mcx-accent); }
.mcx-page { max-width:1680px; margin:0 auto; padding:22px 26px 48px; }
.mcx-head h2 { font-size:1.3rem; font-weight:650; margin:0; letter-spacing:-.01em; }
.mcx-head p { color:var(--mcx-muted); margin:.2rem 0 1.1rem; font-size:.9rem; }
.card { border:1px solid var(--mcx-line); border-radius:12px; box-shadow:0 1px 2px rgba(16,24,40,.04); }
.card-header { background:#fff; font-weight:600; font-size:.9rem; border-bottom:1px solid var(--mcx-line); padding:.7rem 1rem; }
.card-footer { background:#FAFBFC; border-top:1px solid var(--mcx-line); padding:.55rem 1rem; }
.bslib-sidebar-layout { border:1px solid var(--mcx-line) !important; border-radius:12px !important; background:#fff; overflow:visible; }
.bslib-sidebar-layout > .sidebar { background:#fff !important; border-right:1px solid var(--mcx-line) !important; }
.bslib-sidebar-layout > .main { background:var(--mcx-bg); padding:16px; }
.sidebar-content { padding-top:18px !important; }
.mcx-section { margin-bottom:20px; }
.mcx-section-title { font-size:.7rem; font-weight:700; text-transform:uppercase; letter-spacing:.07em; color:#8A90A0; margin-bottom:8px; }
.form-label, .control-label { font-size:.85rem; font-weight:500; color:#374151; margin-bottom:.3rem; }
.form-control, .form-select, .selectize-input { font-size:.88rem; border-color:#D9DCE2; border-radius:8px; }
.form-group, .shiny-input-container { margin-bottom:.8rem; }
.mcx-hint { font-size:.77rem; color:var(--mcx-muted); margin:-4px 0 10px; line-height:1.35; }
.mcx-tip { color:#A0A6B1; margin-left:5px; cursor:help; font-size:.85em; }
.mcx-disabled { opacity:.42; pointer-events:none; filter:grayscale(1); }
.btn { border-radius:8px; font-weight:500; }
.btn-primary { background:var(--mcx-accent); border-color:var(--mcx-accent); }
.btn-outline-primary { color:var(--mcx-accent); border-color:var(--mcx-accent); }
.btn-light { background:#F3F4F6; border-color:#E5E7EB; }
.dl-bar { display:flex; align-items:center; gap:8px; flex-wrap:wrap; font-size:.82rem; color:var(--mcx-muted); }
.dl-bar .dl-title { font-weight:600; color:#4B5563; margin-right:4px; }
.dl-bar .form-group, .dl-bar .shiny-input-container { margin:0 !important; }
.dl-bar .form-control, .dl-bar .form-select, .dl-bar .selectize-input { height:30px; min-height:30px; padding:2px 8px; font-size:.82rem; }
.mcx-stats { display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:12px; margin-bottom:16px; }
.mcx-stat { background:#fff; border:1px solid var(--mcx-line); border-radius:12px; padding:14px 18px; }
.mcx-stat .v { font-size:1.55rem; font-weight:700; letter-spacing:-.02em; }
.mcx-stat .l { font-size:.72rem; color:var(--mcx-muted); text-transform:uppercase; letter-spacing:.06em; font-weight:600; }
.mcx-chip { display:inline-flex; align-items:center; gap:6px; font-size:.78rem; padding:3px 10px; border-radius:999px; background:var(--mcx-accent-soft); color:var(--mcx-accent); font-weight:500; margin-left:6px; white-space:nowrap; }
.mcx-chip.warn { background:#FFF7ED; color:#C2410C; } .mcx-chip.ok { background:#ECFDF5; color:#047857; } .mcx-chip.muted { background:#F3F4F6; color:#6B7280; }
.mcx-empty { text-align:center; padding:70px 20px; color:var(--mcx-muted); } .mcx-empty .i { font-size:2.2rem; color:#C7CBD4; margin-bottom:10px; } .mcx-empty h5 { color:#374151; font-weight:600; }
.mcx-files { padding:0; margin:8px 0 14px; } .mcx-files li { list-style:none; padding:5px 0; font-size:.87rem; display:flex; gap:8px; align-items:baseline; border-bottom:1px dashed #EEF0F3; }
.mcx-files .ok { color:#059669; } .mcx-files .no { color:#DC2626; } .mcx-files .opt { color:#9CA3AF; } .mcx-files code { font-size:.78rem; color:#6B7280; background:none; }
.mcx-needs h6 { font-size:.8rem; font-weight:700; text-transform:uppercase; letter-spacing:.05em; color:#8A90A0; margin-top:6px; }
.mcx-needs ul { padding-left:1.1rem; font-size:.87rem; } .mcx-needs li { margin-bottom:6px; }
.mcx-seg .shiny-options-group { display:inline-flex; flex-wrap:nowrap; border:1px solid #D9DCE2; border-radius:8px; overflow:hidden; background:#fff; }
.mcx-seg .radio-inline, .mcx-seg .form-check-inline { margin:0 !important; padding:0 !important; }
.mcx-seg input[type=radio] { position:absolute; opacity:0; pointer-events:none; }
.mcx-seg label span { display:block; padding:5px 9px; font-size:.81rem; white-space:nowrap; border-right:1px solid #E5E7EB; cursor:pointer; color:#374151; }
.mcx-seg .radio-inline:last-child label span, .mcx-seg label:last-child span { border-right:none; }
.mcx-seg input:checked + span { background:var(--mcx-accent); color:#fff; }
.mcx-seg input:disabled + span { opacity:.38; cursor:not-allowed; }
.mcx-selinfo { background:var(--mcx-accent-soft); color:#1E3A8A; border-radius:8px; padding:8px 12px; font-size:.86rem; margin:10px 0; }
.mcx-colrow { display:flex; align-items:center; gap:8px; margin-bottom:5px; font-size:.84rem; }
.mcx-colrow input[type=color] { width:28px; height:22px; border:1px solid #D9DCE2; border-radius:5px; padding:0; background:none; cursor:pointer; }
.swatch { display:inline-block; width:10px; height:10px; border-radius:3px; margin-right:7px; vertical-align:middle; }
.accordion { --bs-accordion-border-color:var(--mcx-line); --bs-accordion-btn-padding-y:.55rem; }
.accordion-button { font-size:.86rem; font-weight:600; color:#374151; }
.accordion-button:not(.collapsed) { background:#F8F9FB; color:var(--mcx-ink); box-shadow:none; }
table.dataTable { font-size:.85rem; } table.dataTable thead th { font-weight:600; color:#374151; border-bottom:1px solid var(--mcx-line) !important; }
table.dataTable tbody tr.selected > * { box-shadow: inset 0 0 0 9999px rgba(53,87,214,.14) !important; color:var(--mcx-ink) !important; }
pre.shiny-text-output { background:#F8F9FB; border:1px solid var(--mcx-line); border-radius:8px; font-size:.8rem; white-space:pre-wrap; }
.mcx-genes { font-size:.8rem; color:#4B5563; background:#F8F9FB; border-radius:8px; padding:8px 10px; max-height:120px; overflow:auto; margin-bottom:12px; }
.mcx-caveat { font-size:.78rem; color:var(--mcx-muted); padding:6px 2px 0; }
.nav-underline .nav-link { font-size:.88rem; }
html.shiny-busy::before { content:''; position:fixed; top:0; left:0; height:3px; width:40%; background:linear-gradient(90deg,transparent,var(--mcx-accent),transparent); animation:mcxbusy 1s linear infinite; z-index:10000; }
@keyframes mcxbusy { from { transform:translateX(-100%); } to { transform:translateX(260%); } }
.shiny-notification { border-radius:10px; border:1px solid var(--mcx-line); box-shadow:0 8px 24px rgba(16,24,40,.12); font-size:.88rem; }
"

JS <- "
(function wait() {
  if (!(window.Shiny && Shiny.inputBindings && Shiny.addCustomMessageHandler)) { setTimeout(wait, 30); return; }
  var b = new Shiny.InputBinding();
  $.extend(b, {
    find: function(scope) { return $(scope).find('input.mcx-color'); },
    getValue: function(el) { return el.value; },
    setValue: function(el, v) { el.value = v; },
    subscribe: function(el, cb) { $(el).on('change.mcxColor', function() { cb(false); }); },
    unsubscribe: function(el) { $(el).off('.mcxColor'); }
  });
  Shiny.inputBindings.register(b, 'mcx.color');
  Shiny.addCustomMessageHandler('mcx-disable', function(m) {
    (m.ids || []).forEach(function(id) {
      var el = document.getElementById(id); if (!el) return;
      var grp = $(el).closest('.shiny-input-container, .form-group, .bslib-input-switch');
      if (!grp.length) grp = $(el);
      grp.toggleClass('mcx-disabled', !!m.disable);
      $(el).find('input,select,button,textarea').addBack().prop('disabled', !!m.disable);
      if (el.selectize) { m.disable ? el.selectize.disable() : el.selectize.enable(); }
    });
  });
  Shiny.addCustomMessageHandler('mcx-disable-choice', function(m) {
    $('input[name=\"' + m.id + '\"][value=\"' + m.value + '\"]').prop('disabled', !!m.disable);
  });
})();
"

# ============================================================== UI

ui <- page_navbar(
  title = span(span(class = "mcx-logo"), "MCexplorer"),
  id = "nav", fillable = FALSE, window_title = "MCexplorer",
  theme = bs_theme(version = 5, primary = "#3557D6", secondary = "#6B7280", "body-bg" = "#F4F5F7",
                   base_font = font_collection("Inter", "system-ui", "-apple-system", "Segoe UI", "Roboto", "Helvetica Neue", "Arial", "sans-serif"),
                   code_font = font_collection("SFMono-Regular", "Menlo", "Consolas", "monospace"), font_scale = 0.95),
  header = tags$head(tags$style(HTML(CSS)), tags$script(HTML(JS))),

  # ---------------------------------------------------------------- Data
  nav_panel("Data", value = "data",
    div(class = "mcx-page",
      page_head("Data", "Open the files from the metacell gene-module step. Everything else is prepared automatically."),
      uiOutput("overview_ui"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Open a dataset"),
          navset_underline(id = "src",
            nav_panel("From a folder", value = "folder", br(),
              textInput("folder", "Folder with your metacell files", value = DATA_PATH, width = "100%",
                        placeholder = "/full/path/to/the/folder"),
              actionButton("scan", "Look in folder", icon = icon("magnifying-glass"), class = "btn-outline-primary"),
              uiOutput("scan_ui")),
            nav_panel("Upload files", value = "upload", br(),
              fileInput("up_files", NULL, multiple = TRUE, width = "100%", buttonLabel = "Choose files",
                        placeholder = "Select all files at once (Ctrl / Cmd click)",
                        accept = c(".rd", ".rda", ".RData", ".rds", ".RDS", ".csv", ".tsv", ".txt")),
              uiOutput("up_ui")))),
        card(class = "mcx-needs", card_header("What the app needs"),
          h6("Required (4 files)"),
          tags$ul(
            tags$li(strong("Gene modules"), " (file name contains ", code("genelist"), ")"),
            tags$li(strong("Module scores"), " (", code("metadata"), ")"),
            tags$li(strong("Gene expression"), " (", code("counts"), ")"),
            tags$li(strong("Module heatmap"), " (", code("HMlist"), ")")),
          p(class = "mcx-hint", "These are the .rd files saved by the metacell modules scoring script. When uploading, file names do not matter: the app recognises each file by its content."),
          h6("Optional"),
          tags$ul(
            tags$li(strong("Cell information"), ": a Seurat object (.RDS) or a table of cells (.csv / .tsv) with a column of metacell names. ",
                    "Adds cell types, samples, conditions and cell counts. Large objects are read in the background while you work, and only once."),
            tags$li(strong("Metacell UMAP"), " (.csv with metacell names and two coordinate columns). Without it, the app computes a UMAP from the module scores.")),
          p(class = "mcx-hint", "Uploading a large Seurat object through the browser is slow; a small table of cells (cell, metacell, sample, condition...) works just as well.")))
    )),

  # ---------------------------------------------------------------- Module heatmap (saved)
  nav_panel("Module heatmap", value = "heat1",
    tab_body("Module heatmap", "The gene-module heatmap as saved by the pipeline. Columns are metacells grouped into blocks; rows are gene modules.",
      tagList(
        section("Display",
          input_switch("h1_blocknames", lab("Number blocks left to right", "Blocks are numbered by their position in the heatmap (B1 = leftmost). These are the numbers used everywhere in the app and in your scripts' 'block' variable."), TRUE),
          selectizeInput("h1_bars", "Label bars above the heatmap", choices = NULL, multiple = TRUE,
                         options = list(placeholder = "None", plugins = list("remove_button")))),
        section("Size", sliderInput("h1_height", "Height on screen (px)", 300, 2400, 700, step = 50))),
      plot_card("h1_plot", "h1dl", 12, 7))),

  # ---------------------------------------------------------------- Heatmap + genes
  nav_panel("Heatmap + genes", value = "heat2",
    tab_body("Heatmap + genes", "Module scores per metacell with gene expression underneath. Choose which genes and blocks to show.",
      tagList(
        section("Genes underneath",
          radioButtons("g_src", NULL, c("Top genes of each block" = "top", "Most variable genes" = "var",
                                        "My own list" = "custom", "Original app list (T/NK cells)" = "orig")),
          conditionalPanel("input.g_src == 'top'", numericInput("g_top_n", "Genes per block", 2, min = 1, max = 10)),
          conditionalPanel("input.g_src == 'var'", numericInput("g_var_n", "Number of genes", 40, min = 5, max = 300)),
          conditionalPanel("input.g_src == 'custom'", textAreaInput("g_custom", NULL, rows = 4, placeholder = "CD3D, CD14, MS4A1 ...")),
          uiOutput("g_status")),
        section("Blocks",
          radioButtons("blk_mode", NULL, c("All, heatmap order" = "all", "Grouped by cell type" = "type", "Choose blocks" = "pick")),
          conditionalPanel("input.blk_mode == 'type'", selectInput("blk_type_col", NULL, choices = NULL)),
          conditionalPanel("input.blk_mode == 'pick'",
            selectizeInput("blk_pick", NULL, choices = NULL, multiple = TRUE, options = list(plugins = list("remove_button"), placeholder = "Pick blocks")),
            hint("Shown in the order you pick them.")),
          hint("Also used by 'Module by block', 'Compare' and 'Composition'.")),
        section("Label bar", selectInput("h2_bar", NULL, choices = c("None" = "none"))),
        actionButton("h2_go", "Update heatmap", icon = icon("rotate"), class = "btn-primary w-100"),
        br(), br(),
        accordion(open = FALSE,
          accordion_panel("More options", icon = icon("sliders"),
            input_switch("h2_genez", lab("Scale each gene", "Colour each gene relative to its own range (z-score) instead of absolute expression."), FALSE),
            input_switch("h2_sep", "Lines between blocks", TRUE),
            sliderInput("h2_gprop", "Space for genes", 0, 0.8, 0.4, step = 0.05),
            sliderInput("h2_lab", "Row label size", 3, 14, 7, step = 0.5),
            sliderInput("h2_height", "Height on screen (px)", 400, 3600, 1000, step = 50),
            input_switch("h2_hover", lab("Hover to see values", "Interactive version. Slower with many metacells."), FALSE),
            fileInput("mc_order_file", lab("Custom metacell order (.csv)", "First column lists metacell names in the order to plot. Metacells not listed are left out."), accept = c(".csv", ".txt")),
            input_switch("h2_custom_order", "Use custom metacell order", FALSE)))),
      card(full_screen = TRUE, card_body(uiOutput("h2_out")), card_footer(dl_bar("h2dl", 12, 10))))),

  # ---------------------------------------------------------------- Module by block
  nav_panel("Module by block", value = "byblock",
    tab_body("Module by block", "How strongly one gene module is active in each block. Each dot is a metacell.",
      tagList(
        section("Gene module", selectInput("bb_mod", NULL, choices = NULL), uiOutput("bb_genes")),
        section("Display",
          selectInput("bb_fill", lab("Colour blocks by", "Each block is coloured by the label most of its metacells carry."), choices = c("Nothing" = "none")),
          sliderInput("bb_height", "Height on screen (px)", 300, 1600, 520, step = 20))),
      plot_card("bb_plot", "bbdl", 11, 5))),

  # ---------------------------------------------------------------- Annotate
  nav_panel("Annotate", value = "annot",
    div(class = "mcx-page",
      page_head("Annotate", "Pick metacells, give them a name, repeat. Your names can then be compared, plotted and exported."),
      needs_data(
        layout_columns(col_widths = c(8, 4),
          card(card_header(span("1  Choose metacells")),
            card_body(
              seg(radioButtons("sel_mode", NULL, inline = TRUE, c("By block" = "block", "By cell type" = "type", "By range" = "range", "By name" = "names"))),
              conditionalPanel("input.sel_mode == 'block'",
                hint("Click rows to select blocks. Click again to deselect."),
                DTOutput("sel_blk_tbl"),
                div(class = "mt-2", actionButton("sel_all", "Select all", class = "btn-light btn-sm"), actionButton("sel_none", "Clear", class = "btn-light btn-sm"))),
              conditionalPanel("input.sel_mode == 'type'",
                layout_columns(col_widths = c(5, 7),
                  selectInput("sel_type_col", "Label", choices = NULL),
                  selectizeInput("sel_type_val", "Values", choices = NULL, multiple = TRUE, options = list(plugins = list("remove_button"))))),
              conditionalPanel("input.sel_mode == 'range'",
                hint("Metacells between two names, in left-to-right heatmap order."),
                layout_columns(col_widths = c(6, 6),
                  selectizeInput("sel_from", "From", choices = NULL), selectizeInput("sel_to", "To", choices = NULL))),
              conditionalPanel("input.sel_mode == 'names'",
                textAreaInput("sel_names", NULL, rows = 4, width = "100%", placeholder = "Paste metacell names separated by commas, spaces or new lines")),
              accordion(open = FALSE,
                accordion_panel("Leave some out (optional)", icon = icon("filter"),
                  hint("Remove a range of metacells (heatmap order) from the selection above."),
                  layout_columns(col_widths = c(6, 6),
                    selectizeInput("ex_from", "From", choices = NULL), selectizeInput("ex_to", "To", choices = NULL)))),
              uiOutput("sel_info"))),
          card(card_header("2  Name them"),
            card_body(
              textInput("ann_name", NULL, placeholder = "e.g. Inflammatory macrophages", width = "100%"),
              div(class = "d-flex gap-2", actionButton("ann_assign", "Assign name", icon = icon("tag"), class = "btn-primary flex-grow-1"),
                  actionButton("ann_undo", NULL, icon = icon("rotate-left"), class = "btn-light", title = "Undo last change")),
              hr(),
              div(class = "mcx-section-title", "Your annotation"),
              DTOutput("ann_tbl"),
              hr(),
              div(class = "d-flex flex-wrap gap-2",
                downloadButton("ann_dl_csv", "Annotation (.csv)", class = "btn-sm btn-outline-primary"),
                downloadButton("ann_dl_r", "Relabel code (.R)", class = "btn-sm btn-outline-primary"),
                downloadButton("ann_dl_blk", "Block summary (.csv)", class = "btn-sm btn-outline-primary")),
              br(),
              accordion(open = FALSE,
                accordion_panel("Other options", icon = icon("ellipsis"),
                  selectInput("ann_init_col", lab("Start from an existing labelling", "Replaces your annotation with a cell-type label from the cell information file. Undo is available."), choices = c("None" = "none")),
                  actionButton("ann_init", "Use this labelling", class = "btn-light btn-sm"),
                  hr(),
                  fileInput("ann_file", "Load a saved annotation (.csv)", accept = ".csv"),
                  uiOutput("ann_load_msg"),
                  hr(),
                  input_switch("autosave", "Save every change automatically", TRUE),
                  uiOutput("autosave_ui")))))),
        accordion(open = FALSE, accordion_panel("History of changes", icon = icon("clock-rotate-left"), DTOutput("ann_log"))))
    )),

  # ---------------------------------------------------------------- Compare
  nav_panel("Compare", value = "compare",
    tab_body("Compare groups", "Compare a gene module, a gene, a set of genes or a metacell metric between groups of metacells.",
      tagList(
        section("Groups",
          selectInput("cg_by", "Compare", choices = NULL),
          selectizeInput("cg_groups", lab("Groups to show", "Drag-free ordering: groups appear in the order listed here. Remove a group to leave it out."),
                         choices = NULL, multiple = TRUE, options = list(plugins = list("remove_button")))),
        section("Measure",
          seg(radioButtons("cg_what", NULL, inline = TRUE, c("Module" = "module", "Gene" = "gene", "Gene set" = "set", "Metric" = "metric"))),
          conditionalPanel("input.cg_what == 'module'", selectInput("cg_mod", NULL, choices = NULL)),
          conditionalPanel("input.cg_what == 'gene'", selectizeInput("cg_gene", NULL, choices = NULL)),
          conditionalPanel("input.cg_what == 'set'",
            textAreaInput("cg_set", NULL, rows = 3, placeholder = "Two or more genes"),
            hint("Combined score of these genes, scaled across the metacells shown.")),
          conditionalPanel("input.cg_what == 'metric'", selectInput("cg_metric", NULL, choices = NULL))),
        section("Statistics",
          selectInput("cg_stats", NULL, c("No statistics" = "none", "Show p-values" = "p", "Show stars" = "stars")),
          radioButtons("cg_cmp", "Comparisons", c("All pairs" = "all", "Each group against one" = "ref", "Pick pairs" = "pick")),
          selectInput("cg_ref", "Reference group", choices = NULL),
          selectizeInput("cg_pairs", "Pairs", choices = NULL, multiple = TRUE, options = list(plugins = list("remove_button"))),
          selectInput("cg_adj", "Correction for multiple tests", c("Bonferroni" = "bonferroni", "Benjamini-Hochberg (FDR)" = "BH", "Holm" = "holm", "None" = "none")),
          input_switch("cg_hidens", "Hide non-significant", FALSE)),
        accordion(open = FALSE,
          accordion_panel("Colours", icon = icon("palette"), uiOutput("cg_colors"), actionButton("cg_colreset", "Reset colours", class = "btn-light btn-sm")),
          accordion_panel("Appearance", icon = icon("sliders"),
            input_switch("cg_legend", "Legend", TRUE),
            input_switch("cg_mean", lab("Overall mean line", "Dashed line at the mean of all metacells shown."), TRUE),
            sliderInput("cg_pt", "Dot size", 0, 3, 0.8, step = 0.1),
            sliderInput("cg_txt", "Text size", 8, 24, 13, step = 1),
            sliderInput("cg_statsize", "Statistics label size", 2, 8, 3.5, step = 0.25),
            sliderInput("cg_height", "Height on screen (px)", 300, 1600, 560, step = 20)))),
      plot_card("cg_plot", "cgdl", 7, 6),
      card(card_header(div(class = "d-flex justify-content-between align-items-center", "Statistics",
                           downloadButton("cg_dl_stats", "Table (.csv)", class = "btn-sm btn-light"))),
           DTOutput("cg_stats_tbl"),
           div(class = "mcx-caveat", icon("circle-info"), " Each dot is a metacell. p-values treat metacells as independent, so use them to rank differences, not as sample-level evidence. Kruskal-Wallis across all groups, pairwise Wilcoxon tests for the brackets.")))),

  # ---------------------------------------------------------------- Composition
  nav_panel("Composition", value = "comp",
    div(class = "mcx-page",
      page_head("Composition", "Which samples or conditions the cells of each block or group come from. Counts are cells, not metacells."),
      needs_data(
        conditionalPanel("!output.has_comp", empty_state("table-cells", "No cell information yet",
          "Composition needs a cell information file with sample or condition columns. Add it on the Data page (or wait if it is still loading).")),
        conditionalPanel("output.has_comp",
          layout_sidebar(border = TRUE, fillable = FALSE,
            sidebar = sidebar(width = 320,
              section("Show", selectInput("cp_col", NULL, choices = NULL)),
              section("For each", selectInput("cp_by", NULL, choices = NULL)),
              section("Scale",
                radioButtons("cp_scale", NULL, c("% of cells" = "raw", "Adjusted for group size" = "norm")),
                hint("Adjusted: each sample / condition is first divided by its total number of cells, so groups with more cells sequenced do not dominate.")),
              accordion(open = FALSE,
                accordion_panel("Colours", icon = icon("palette"), uiOutput("cp_colors")),
                accordion_panel("Appearance", icon = icon("sliders"),
                  input_switch("cp_flip", "Horizontal bars", FALSE),
                  sliderInput("cp_height", "Height on screen (px)", 300, 1600, 520, step = 20)))),
            plot_card("cp_plot", "cpdl", 10, 5)))))),

  # ---------------------------------------------------------------- Find genes
  nav_panel("Find genes", value = "genes",
    div(class = "mcx-page",
      page_head("Find genes", "Look up genes and modules, and find the genes that set blocks apart."),
      needs_data(
        navset_underline(id = "fg_tabs",
          nav_panel("A gene", br(),
            layout_columns(col_widths = c(4, 8),
              card(card_body(
                selectizeInput("fg_gene", "Gene", choices = NULL),
                uiOutput("fg_info"),
                actionButton("fg_add", "Show on heatmap", icon = icon("plus"), class = "btn-light btn-sm"))),
              plot_card("fg_plot", "fgdl", 10, 4))),
          nav_panel("Top genes of blocks", br(),
            layout_columns(col_widths = c(4, 8),
              card(card_body(
                selectizeInput("tg_a", "Blocks", choices = NULL, multiple = TRUE, options = list(plugins = list("remove_button"))),
                radioButtons("tg_vs", "Compared with", c("All other metacells" = "rest", "Chosen blocks" = "blocks")),
                selectizeInput("tg_b", NULL, choices = NULL, multiple = TRUE, options = list(plugins = list("remove_button"), placeholder = "Blocks to compare with")),
                numericInput("tg_n", "Number of genes", 50, min = 5, max = 1000),
                div(class = "d-flex flex-wrap gap-2",
                  actionButton("tg_add", "Show top 20 on heatmap", icon = icon("plus"), class = "btn-light btn-sm"),
                  downloadButton("tg_dl", "Table (.csv)", class = "btn-light btn-sm")),
                hint("Ranked by log2 fold change of mean log2(fraction + 1e-5), the usual metacell convention."))),
              card(DTOutput("tg_tbl")))),
          nav_panel("Top modules of a block", br(),
            layout_columns(col_widths = c(4, 8),
              card(card_body(
                selectInput("tm_block", "Block", choices = NULL),
                sliderInput("tm_thr", lab("Minimum average module score", "Average z-score of the module across the block's metacells. The original app used 0.8."), 0, 3, 0.8, step = 0.1))),
              card(DTOutput("tm_tbl")))),
          nav_panel("Genes varying within blocks", br(),
            layout_columns(col_widths = c(4, 8),
              card(card_body(
                selectizeInput("vg_blocks", "Blocks", choices = NULL, multiple = TRUE, options = list(plugins = list("remove_button"))),
                actionButton("vg_go", "Find genes", icon = icon("magnifying-glass"), class = "btn-primary"),
                uiOutput("vg_note"),
                hint("Chi-square test per gene: are its UMIs spread across these metacells in proportion to metacell size? Top genes vary most between the metacells of the chosen blocks."),
                div(class = "d-flex flex-wrap gap-2",
                  actionButton("vg_add", "Show top 30 on heatmap", icon = icon("plus"), class = "btn-light btn-sm"),
                  downloadButton("vg_dl", "Table (.csv)", class = "btn-light btn-sm")))),
              card(DTOutput("vg_tbl")))),
          nav_panel("All modules", br(),
            card(card_header(div(class = "d-flex justify-content-between align-items-center", "Gene modules",
                                 downloadButton("mod_dl", "Table (.csv)", class = "btn-sm btn-light"))), DTOutput("mod_tbl"))),
          nav_panel("All blocks", br(),
            card(card_header(div(class = "d-flex justify-content-between align-items-center", "Blocks",
                                 downloadButton("blk_dl", "Table (.csv)", class = "btn-sm btn-light"))), DTOutput("blk_tbl"))))))),

  # ---------------------------------------------------------------- UMAP
  nav_panel("UMAP", value = "umap",
    tab_body("UMAP", "Metacells in two dimensions, coloured by anything in the app.",
      tagList(
        section("Colour by",
          selectInput("um_by", NULL, choices = NULL),
          conditionalPanel("input.um_by == 'module'", selectInput("um_mod", NULL, choices = NULL)),
          conditionalPanel("input.um_by == 'gene'", selectizeInput("um_gene", NULL, choices = NULL))),
        section("Display",
          input_switch("um_labels", "Group labels", TRUE),
          sliderInput("um_pt", "Dot size", 0.2, 5, 1.4, step = 0.1),
          sliderInput("um_height", "Height on screen (px)", 300, 1600, 640, step = 20)),
        accordion(open = FALSE,
          accordion_panel("UMAP settings", icon = icon("gear"),
            uiOutput("um_source"),
            sliderInput("um_nn", "Neighbours", 5, 100, 15, step = 1),
            sliderInput("um_md", "Minimum distance", 0, 1, 0.3, step = 0.05),
            actionButton("um_recompute", "Recompute", icon = icon("rotate"), class = "btn-light btn-sm")))),
      plot_card("um_plot", "umdl", 8, 7))),

  # ---------------------------------------------------------------- Save
  nav_panel("Save", value = "save",
    div(class = "mcx-page",
      page_head("Save", "Keep your work: one zip with everything, or a session file to pick up exactly where you left off."),
      needs_data(
        layout_columns(col_widths = c(4, 4, 4),
          card(card_header("Export everything"), card_body(
            p(class = "mcx-hint", "Annotation (.csv + relabel code + history), all tables, every figure you have looked at in this session (PDF, at the sizes set under each figure), statistics, and a session file."),
            downloadButton("dl_all", "Download zip", class = "btn-primary"))),
          card(card_header("Save session"), card_body(
            p(class = "mcx-hint", "Your annotation, its history and all settings, in one small file."),
            downloadButton("dl_session", "Download session", class = "btn-outline-primary"))),
          card(card_header("Restore session"), card_body(
            p(class = "mcx-hint", "Open the matching dataset first, then load the session file."),
            fileInput("session_file", NULL, accept = ".rds"), uiOutput("session_msg"))))))),

  nav_spacer(),
  nav_item(uiOutput("nav_status"))
)

# inputs stored in a session file
RESTORE <- c(h1_blocknames = "switch", h1_bars = "multi", h1_height = "slider",
  g_src = "radio", g_top_n = "numeric", g_var_n = "numeric", g_custom = "text", blk_mode = "radio", blk_type_col = "select",
  blk_pick = "multi", h2_bar = "select", h2_genez = "switch", h2_sep = "switch", h2_gprop = "slider", h2_lab = "slider",
  h2_height = "slider", h2_hover = "switch", bb_mod = "select", bb_fill = "select", bb_height = "slider", autosave = "switch",
  cg_by = "select", cg_groups = "multi", cg_what = "radio", cg_mod = "select", cg_set = "text", cg_metric = "select",
  cg_stats = "select", cg_cmp = "radio", cg_ref = "select", cg_pairs = "multi", cg_adj = "select", cg_hidens = "switch",
  cg_legend = "switch", cg_mean = "switch", cg_pt = "slider", cg_txt = "slider", cg_statsize = "slider", cg_height = "slider",
  cp_col = "select", cp_by = "select", cp_scale = "radio", cp_flip = "switch", cp_height = "slider",
  um_by = "select", um_mod = "select", um_labels = "switch", um_pt = "slider", um_height = "slider")

# ============================================================== server

server <- function(input, output, session) {

  D <- reactiveVal(NULL)       # core dataset
  CELLS <- reactiveVal(NULL)   # per-metacell summary of the cell information
  UMAP <- reactiveVal(NULL)    # list(xy, source)
  rv <- reactiveValues(scan = NULL, up = NULL, ann = NULL, hist = list(), log = NULL, colors = character(0),
                       job = NULL, cell_status = list(state = "none"), vg = NULL, session_msg = NULL, ann_msg = NULL,
                       umap_msg = NULL, pending_umap = FALSE, nonce = 0)
  colmem <- new.env()          # automatic colour per name, stable for the session
  cache_env <- new.env()       # cache contents of the open dataset (not reactive)
  figs <- new.env(); tables <- new.env()

  disable <- function(ids, state = TRUE) session$sendCustomMessage("mcx-disable", list(ids = as.list(ids), disable = state))
  disable_choice <- function(id, value, state = TRUE) session$sendCustomMessage("mcx-disable-choice", list(id = id, value = value, disable = state))

  output$has_data <- reactive(!is.null(D()))
  output$has_comp <- reactive(length(CELLS()$comp) > 0)
  outputOptions(output, "has_data", suspendWhenHidden = FALSE)
  outputOptions(output, "has_comp", suspendWhenHidden = FALSE)

  # ------------------------------------------------------------ colours
  pal_for <- function(lv) {
    lv <- as.character(lv); user <- rv$colors
    vapply(lv, function(x) {
      if (!is.na(user[x])) return(unname(user[x]))
      if (is.na(x) || x %in% NO_LABEL) return(GREY)
      if (!exists(x, envir = colmem, inherits = FALSE)) {
        n <- length(ls(colmem, all.names = TRUE))
        assign(x, if (n < length(PAL)) PAL[n + 1] else grDevices::hcl.colors(n + 1, "Dynamic")[n + 1], envir = colmem)
      }
      get(x, envir = colmem)
    }, "", USE.NAMES = TRUE)
  }

  # ------------------------------------------------------------ dataset values
  label_cols <- reactive(names(CELLS()$labels) %||% character(0))
  comp_cols  <- reactive(names(CELLS()$comp) %||% character(0))
  primary_label <- reactive({
    lc <- label_cols(); if (!length(lc)) return(NULL)
    p <- lc[grepl("label|type|annot|lineage|state", lc, ignore.case = TRUE) & !grepl("block|module", lc, ignore.case = TRUE)]
    (p %||% lc)[1]
  })

  mv <- function(key, mcs = NULL) {
    b <- D(); cl <- CELLS(); if (is.null(mcs)) mcs <- b$mcs
    if (key == "ann") return(rv$ann$annotation[match(mcs, rv$ann$metacell)])
    if (key == "block") return(paste0("B", b$block[mcs]))
    type <- sub(":.*$", "", key); col <- sub("^[^:]*:", "", key); i <- match(mcs, b$mcs)
    switch(type,
      lab = cl$labels[[col]][i], maj = cl$comp_major[[col]][i], pur = cl$comp_purity[[col]][i],
      num = switch(col, n_genes = b$nfeat[i], n_cells = cl$n_cells[i], total_umis = cl$total_umis[i], cl$nums[[col]][i]),
      NULL)
  }
  key_label <- function(key) {
    if (key == "ann") return("Your annotation"); if (key == "block") return("Block")
    col <- sub("^[^:]*:", "", key)
    switch(sub(":.*$", "", key), lab = pretty_name(col), maj = paste("Main", pretty_name(col)), pur = paste("Share of main", pretty_name(col)),
           num = switch(col, n_genes = "Genes detected", n_cells = "Cells", total_umis = "UMIs", pretty_name(col)), col)
  }
  level_order <- function(key, v) {
    v <- unique(as.character(v[!is.na(v)]))
    if (key == "block") return(intersect(paste0("B", blocks_shown()), v))
    col <- sub("^[^:]*:", "", key); lv <- CELLS()$levels[[col]]
    if (!is.null(lv)) return(c(intersect(lv, v), setdiff(v, lv)))
    mixed_sort(v)
  }
  group_choices <- function(annotation = TRUE, block = TRUE) {
    cl <- CELLS(); out <- list()
    if (annotation) out[["Your work"]] <- c("Your annotation" = "ann")
    if (block) out[["Heatmap"]] <- c("Block" = "block")
    if (length(cl$labels)) out[["Cell type labels"]] <- setNames(paste0("lab:", names(cl$labels)), pretty_name(names(cl$labels)))
    if (length(cl$comp)) out[["Most common per metacell"]] <- setNames(paste0("maj:", names(cl$comp)), paste("Main", pretty_name(names(cl$comp))))
    out
  }
  metric_choices <- function() {
    cl <- CELLS(); out <- c("Genes detected" = "num:n_genes")
    if (!is.null(cl$n_cells)) out <- c(out, "Cells" = "num:n_cells")
    if (!is.null(cl$total_umis)) out <- c(out, "UMIs" = "num:total_umis")
    if (length(cl$nums)) out <- c(out, setNames(paste0("num:", names(cl$nums)), paste("Mean", pretty_name(names(cl$nums)))))
    if (length(cl$comp_purity)) out <- c(out, setNames(paste0("pur:", names(cl$comp_purity)), paste("Share of main", pretty_name(names(cl$comp_purity)))))
    out
  }
  block_major <- function(key) {
    b <- D(); v <- mv(key)
    vapply(seq_len(b$nb), function(k) { t <- table(v[b$block == k]); if (length(t)) names(t)[which.max(t)] else NA_character_ }, "")
  }
  block_choices <- function() {
    b <- D(); pl <- primary_label()
    nm <- if (!is.null(pl)) paste0("Block ", seq_len(b$nb), "  (", block_major(paste0("lab:", pl)), ")") else paste("Block", seq_len(b$nb))
    setNames(paste0("B", seq_len(b$nb)), nm)
  }
  mods_choices <- function() { b <- D(); o <- order(b$mods$num); setNames(b$mods$id[o], paste0(b$mods$name[o], "  (", b$mods$n_genes, " genes)")) }

  blocks_shown <- reactive({
    b <- req(D()); mode <- input$blk_mode %||% "all"
    if (mode == "pick") { v <- as.integer(sub("^B", "", input$blk_pick %||% character(0))); v <- v[!is.na(v)]; if (length(v)) return(v) }
    if (mode == "type" && nzchar(input$blk_type_col %||% "")) {
      lab <- block_major(input$blk_type_col); return(order(is.na(lab), lab, seq_len(b$nb)))
    }
    seq_len(b$nb)
  })
  mcs_shown <- reactive({ b <- req(D()); b$mcs[unlist(b$col_ord[blocks_shown()])] })

  # ============================================================ DATA page
  do_scan <- function(dir) {
    dir <- sub("/+$", "", trimws(dir))
    if (!nzchar(dir)) return()
    if (!dir.exists(dir)) { rv$scan <- list(error = paste0("Folder not found: ", dir)); return() }
    s <- scan_folder(dir); rv$scan <- s
    if (length(s$sets) == 1 && identical(dir, sub("/+$", "", DATA_PATH)) && is.null(D())) open_folder_set(names(s$sets)[1])
  }
  observeEvent(input$scan, do_scan(input$folder))
  observeEvent(TRUE, { if (nzchar(DATA_PATH)) do_scan(DATA_PATH) }, once = TRUE)

  file_li <- function(ok, label, path = NULL, optional = FALSE)
    tags$li(span(class = if (ok) "ok" else if (optional) "opt" else "no", icon(if (ok) "circle-check" else if (optional) "circle-minus" else "circle-xmark")),
            span(strong(label), if (!is.null(path)) tagList(br(), code(basename(path)))))

  output$scan_ui <- renderUI({
    s <- rv$scan; if (is.null(s)) return(NULL)
    if (!is.null(s$error)) return(div(class = "alert alert-warning mt-3", s$error))
    if (!length(s$sets)) {
      miss <- if (length(s$incomplete)) paste0(" Found some files but missing: ",
        paste(ROLE_NAMES[setdiff(names(ROLE_KEYS), names(s$incomplete[[1]]))], collapse = ", "), ".") else ""
      return(div(class = "alert alert-warning mt-3", "No complete set of the 4 required files in this folder.", miss))
    }
    stems <- names(s$sets)
    tagList(br(),
      if (length(stems) > 1) radioButtons("set", lab("Dataset", "This folder holds more than one complete set of files."), stems, selected = stems[1]) else NULL,
      uiOutput("scan_files"),
      selectInput("cell_file", lab("Cell information (optional)", "Seurat object or table of cells. Gives cell types, samples and conditions."),
                  c("None" = "", setNames(s$cell_cand, basename(s$cell_cand))), selected = guess_cell_file(s$cell_cand, stems[1]), width = "100%"),
      selectInput("umap_file", "Metacell UMAP", c("Compute automatically" = "", setNames(s$umap_cand, basename(s$umap_cand))), width = "100%"),
      actionButton("open_folder", "Open dataset", icon = icon("arrow-right"), class = "btn-primary"))
  })
  output$scan_files <- renderUI({
    s <- rv$scan; req(length(s$sets)); st <- input$set %||% names(s$sets)[1]; f <- s$sets[[st]]; req(f)
    tags$ul(class = "mcx-files", lapply(names(ROLE_KEYS), function(r) file_li(!is.na(f[r]), ROLE_NAMES[[r]], f[r])))
  })
  observeEvent(input$set, { s <- rv$scan; updateSelectInput(session, "cell_file", selected = guess_cell_file(s$cell_cand, input$set)) })

  open_folder_set <- function(stem, cell_file = NULL, umap_file = NULL) {
    s <- rv$scan; f <- s$sets[[stem]]
    if (is.null(cell_file)) cell_file <- guess_cell_file(s$cell_cand, stem)
    open_dataset(as.list(f), stem, cache_path(s$dir, stem), cell_file, umap_file %||% "")
  }
  observeEvent(input$open_folder, open_folder_set(input$set %||% names(rv$scan$sets)[1], input$cell_file %||% "", input$umap_file %||% ""))

  # uploads: every file is recognised by its content
  observeEvent(input$up_files, {
    fi <- input$up_files
    out <- rv$up %||% list(objs = list(), names = list(), cell = NULL, umap = NULL, unknown = character(0))   # add to earlier uploads
    out$unknown <- character(0)
    withProgress(message = "Reading uploaded files", value = 0, {
      for (i in seq_len(nrow(fi))) {
        nm <- fi$name[i]; p <- fi$datapath[i]; incProgress(1 / nrow(fi), detail = nm)
        if (grepl("\\.(rd|rda|rdata)$", nm, ignore.case = TRUE)) {
          e <- new.env(); objs <- tryCatch(load(p, envir = e), error = function(err) character(0))
          for (o in objs) { r <- rd_role(get(o, envir = e)); if (!is.na(r)) { out$objs[[r]] <- get(o, envir = e); out$names[[r]] <- nm } }
          if (!length(objs) || all(is.na(vapply(objs, function(o) rd_role(get(o, envir = e)), "")))) out$unknown <- c(out$unknown, nm)
        } else if (grepl("umap", nm, ignore.case = TRUE)) out$umap <- list(path = p, name = nm)
        else out$cell <- list(path = p, name = nm)
      }
    })
    rv$up <- out
  })
  output$up_ui <- renderUI({
    u <- rv$up
    if (is.null(u)) return(hint("Select the 4 required files (and optionally a cell information file and a UMAP table) in one go."))
    ok <- all(names(ROLE_KEYS) %in% names(u$objs))
    tagList(
      tags$ul(class = "mcx-files",
        lapply(names(ROLE_KEYS), function(r) file_li(r %in% names(u$objs), ROLE_NAMES[[r]], u$names[[r]])),
        file_li(!is.null(u$cell), "Cell information", u$cell$name, optional = TRUE),
        file_li(!is.null(u$umap), "Metacell UMAP", u$umap$name, optional = TRUE)),
      if (length(u$unknown)) div(class = "alert alert-warning", "Not recognised: ", paste(u$unknown, collapse = ", ")),
      if (!ok) hint("Add the missing files: choose them now, the files above are kept.") else NULL,
      actionButton("open_upload", "Open dataset", icon = icon("arrow-right"), class = if (ok) "btn-primary" else "btn-primary disabled"))
  })
  observeEvent(input$open_upload, {
    u <- rv$up; req(all(names(ROLE_KEYS) %in% names(u$objs)))
    stem <- gsub("_+", "_", gsub("(genelist|metadata|counts|hmlist)", "", tolower(sub("\\.[^.]+$", "", u$names$counts))))
    open_dataset(u$objs, sub("^_|_$", "", stem), NULL, u$cell$path %||% "", u$umap$path %||% "")
  })

  open_dataset <- function(src, name, cache_file, cell_file, umap_file) {
    b <- NULL
    withProgress(message = paste("Opening", name), value = 0, {
      b <- tryCatch(prepare_dataset(src, name, cache_file), error = function(e) {
        mcx_log("ERROR: ", conditionMessage(e)); showNotification(paste("Could not open the dataset:", conditionMessage(e)), type = "error", duration = NULL); NULL })
    })
    if (is.null(b)) return(invisible())
    if (!is.null(rv$job$job) && rv$job$job$is_alive()) rv$job$job$kill()
    rv$job <- NULL; rv$cell_status <- list(state = "none"); CELLS(NULL); UMAP(NULL); rv$vg <- NULL
    rm(list = ls(figs), envir = figs); rm(list = ls(tables), envir = tables)
    cache_env$data <- b$cache; b$cache <- NULL
    D(b)
    if (nzchar(umap_file %||% "")) {
      xy <- tryCatch(read_umap(umap_file, b$mcs), error = function(e) { showNotification(paste("UMAP file:", conditionMessage(e)), type = "warning"); NULL })
      if (!is.null(xy)) UMAP(list(xy = xy, source = paste("From file:", basename(umap_file))))
    } else if (!is.null(cache_env$data$umap)) UMAP(cache_env$data$umap)
    if (nzchar(cell_file %||% "")) start_cells(cell_file)
    nav_select("nav", "heat1")
  }

  # ------------------------------------------------------------ cell information (background)
  start_cells <- function(path) {
    b <- D(); key <- paste0("cells:", file_key(path))
    hit <- cache_env$data[[key]]
    if (!is.null(hit)) { CELLS(hit); rv$cell_status <- list(state = "done", file = basename(path)); mcx_log("Cell information loaded from cache"); return() }
    mcx_log("Reading cell information in the background: ", basename(path), " (", fmt_size(file.size(path)), "). The app is usable meanwhile.")
    if (requireNamespace("callr", quietly = TRUE)) {
      f <- summarise_cells; environment(f) <- globalenv()
      out <- tempfile(fileext = ".rds")
      job <- callr::r_bg(function(f, path, mcs, out) { saveRDS(f(path, mcs), out); TRUE },
                         args = list(f = f, path = path, mcs = b$mcs, out = out), supervise = TRUE)
      rv$job <- list(job = job, out = out, t0 = Sys.time(), path = path, key = key, size = file.size(path))
      rv$cell_status <- list(state = "running", file = basename(path), size = file.size(path), t0 = Sys.time())
    } else {
      res <- withProgress(message = "Reading cell information", tryCatch(summarise_cells(path, b$mcs), error = function(e) e))
      finish_cells(res, key, path)
    }
  }
  finish_cells <- function(res, key, path) {
    if (inherits(res, "error")) {
      msg <- sub("^.*?: ", "", conditionMessage(res)); mcx_log("Cell information failed: ", msg)
      rv$cell_status <- list(state = "error", file = basename(path), msg = msg)
      showNotification(paste("Cell information could not be read:", msg), type = "error", duration = NULL); return()
    }
    b <- D(); CELLS(res); cache_env$data[[key]] <- res; cache_write(b$cache_file, b$cache_key, key, res)
    mcx_log(sprintf("Cell information ready: %s cells, labels: %s, cell info: %s", fmt_n(res$n_total),
                    paste(names(res$labels), collapse = ", ") %||% "none", paste(names(res$comp), collapse = ", ") %||% "none"))
    rv$cell_status <- list(state = "done", file = basename(path))
    showNotification(HTML(paste0("<b>Cell information ready.</b> ", fmt_n(res$n_total), " cells; cell types and composition are now available.")), type = "message")
  }
  observe({
    j <- rv$job; req(j)
    if (j$job$is_alive()) {
      invalidateLater(1000)
      rv$cell_status <- list(state = "running", file = basename(j$path), size = j$size, t0 = j$t0)
      return()
    }
    res <- tryCatch({ j$job$get_result(); readRDS(j$out) }, error = function(e) e)
    rv$job <- NULL
    finish_cells(res, j$key, j$path)
  })

  output$nav_status <- renderUI({
    b <- D(); s <- rv$cell_status
    if (is.null(b)) return(div(class = "mcx-navstatus", span(class = "mcx-chip muted", "No dataset open")))
    if (s$state == "running")
      return(div(class = "mcx-navstatus", span(class = "mcx-chip warn", title = paste("Reading", s$file, "in the background"), icon("spinner", class = "fa-spin"),
                 paste0("Reading cell info ", fmt_secs(difftime(Sys.time(), s$t0, units = "secs"))))))
    tip_txt <- switch(s$state, done = "Cell information loaded", error = paste("Cell information failed:", s$msg), "No cell information")
    div(class = "mcx-navstatus", span(class = paste("mcx-chip", if (s$state == "error") "warn" else ""), title = paste0(b$name, " | ", tip_txt),
        icon(switch(s$state, done = "circle-check", error = "triangle-exclamation", "layer-group")), b$name))
  })

  output$overview_ui <- renderUI({
    b <- D(); if (is.null(b)) return(NULL); cl <- CELLS(); s <- rv$cell_status
    stat <- function(v, l) div(class = "mcx-stat", div(class = "v", v), div(class = "l", l))
    tagList(
      div(class = "mcx-stats",
        stat(fmt_n(length(b$mcs)), "Metacells"),
        stat(if (!is.null(cl$n_total) && isTRUE(cl$cell_level)) fmt_n(cl$n_total) else "-", "Cells"),
        stat(ncol(b$Z), "Gene modules"), stat(b$nb, "Blocks"), stat(fmt_n(nrow(b$counts)), "Genes")),
      card(card_header("Current dataset"), card_body(
        tags$ul(class = "mcx-files",
          tags$li(span(class = "ok", icon("circle-check")), span(strong(b$name), br(), code(if (b$src$counts == "(uploaded)") "uploaded files" else dirname(b$src$counts)))),
          tags$li(span(class = if (s$state == "done") "ok" else "opt", icon(if (s$state == "done") "circle-check" else "circle-minus")),
                  span(strong("Cell information: "),
                       switch(s$state, done = paste0(s$file, " | cell types: ", paste(pretty_name(names(cl$labels)), collapse = ", ") %||% "none",
                                                     " | cell info: ", paste(pretty_name(names(cl$comp)), collapse = ", ") %||% "none"),
                              running = paste0(s$file, " is being read in the background (", fmt_size(s$size), "). Everything else works meanwhile."),
                              error = paste0(s$file, " could not be read: ", s$msg), "none"))),
          tags$li(span(class = if (!is.null(UMAP())) "ok" else "opt", icon(if (!is.null(UMAP())) "circle-check" else "circle-minus")),
                  span(strong("UMAP: "), if (!is.null(UMAP())) UMAP()$source else "computed when you open the UMAP page"))))))
  })

  # ============================================================ when a dataset opens
  observeEvent(D(), {
    b <- D()
    rv$ann <- data.frame(metacell = b$mcs, annotation = "undefined", stringsAsFactors = FALSE)
    rv$hist <- list(); rv$log <- data.frame(time = character(0), action = character(0), metacells = integer(0), name = character(0))
    mc <- mods_choices(); bc <- block_choices()
    for (id in c("bb_mod", "cg_mod", "um_mod")) updateSelectInput(session, id, choices = mc)
    updateSelectizeInput(session, "blk_pick", choices = bc, selected = character(0))
    for (id in c("tg_a", "tg_b", "vg_blocks")) updateSelectizeInput(session, id, choices = bc, selected = if (id == "tg_a") "B1" else character(0))
    updateSelectInput(session, "tm_block", choices = bc)
    for (id in c("sel_from", "sel_to", "ex_from", "ex_to")) updateSelectizeInput(session, id, choices = c("", b$mc_order_ht), server = TRUE)
    first_gene <- b$block_top[[1]][1] %||% rownames(b$counts)[1]
    for (id in c("cg_gene", "fg_gene", "um_gene")) updateSelectizeInput(session, id, choices = rownames(b$counts), selected = first_gene, server = TRUE)
    refresh_choices()
    if (file.exists(autosave_file())) rv$ann_msg <- "autosave"
  }, priority = 10)

  observeEvent(CELLS(), refresh_choices(), ignoreNULL = FALSE, ignoreInit = TRUE)

  refresh_choices <- function() {
    b <- D(); req(b); cl <- CELLS(); gc_all <- group_choices(); pl <- primary_label()
    lab_only <- if (length(cl$labels)) setNames(paste0("lab:", names(cl$labels)), pretty_name(names(cl$labels))) else character(0)
    bars <- c(gc_all, list("Metacell" = c("Cells" = "num:n_cells", "Genes detected" = "num:n_genes")[c(!is.null(cl$n_cells), TRUE)]))
    updateSelectizeInput(session, "h1_bars", choices = bars, selected = isolate(input$h1_bars) %||% (if (!is.null(pl)) paste0("lab:", pl)))
    updateSelectInput(session, "h2_bar", choices = c(list("None" = c("None" = "none")), gc_all),
                      selected = if (!is.null(pl)) paste0("lab:", pl) else isolate(input$h2_bar) %||% "none")
    updateSelectInput(session, "blk_type_col", choices = lab_only)
    disable_choice("blk_mode", "type", !length(lab_only))
    updateSelectInput(session, "bb_fill", choices = c("Nothing" = "none", lab_only), selected = if (!is.null(pl)) paste0("lab:", pl) else "none")
    updateSelectInput(session, "sel_type_col", choices = c(lab_only, if (length(cl$comp)) setNames(paste0("maj:", names(cl$comp)), paste("Main", pretty_name(names(cl$comp))))))
    disable_choice("sel_mode", "type", !(length(cl$labels) || length(cl$comp)))
    updateSelectInput(session, "ann_init_col", choices = c("None" = "none", lab_only))
    updateSelectInput(session, "cg_by", choices = gc_all, selected = if (!is.null(pl)) paste0("lab:", pl) else "block")
    updateSelectInput(session, "cg_metric", choices = metric_choices())
    cc <- comp_cols() %||% character(0)
    updateSelectInput(session, "cp_col", choices = if (length(cc)) setNames(cc, pretty_name(cc)) else character(0),
                      selected = if (length(cc)) c(grep("condition|group|disease|treat", cc, value = TRUE, ignore.case = TRUE), cc)[1])
    updateSelectInput(session, "cp_by", choices = group_choices()[c("Heatmap", "Your work", "Cell type labels")[c(TRUE, TRUE, length(cl$labels) > 0)]])
    updateSelectInput(session, "um_by", choices = c(gc_all, list("Expression" = c("Gene module" = "module", "Gene" = "gene")), list("Metacell" = metric_choices())),
                      selected = if (!is.null(pl)) paste0("lab:", pl) else "block")
    updateSelectizeInput(session, "blk_pick", choices = block_choices(), selected = isolate(input$blk_pick))
    for (id in c("tg_a", "tg_b", "vg_blocks")) updateSelectizeInput(session, id, choices = block_choices(),
                                                                     selected = isolate(input[[id]]) %||% (if (id == "tg_a") "B1" else character(0)))
    updateSelectInput(session, "tm_block", choices = block_choices(), selected = isolate(input$tm_block) %||% "B1")
    disable("vg_go", is.null(cl$total_umis))
  }

  # ============================================================ MODULE HEATMAP (saved)
  h1_draw <- reactive({
    b <- req(D()); bars <- input$h1_bars
    bars <- Filter(function(k) !is.null(mv(k, b$mcs[1])), bars %||% character(0))
    if (!length(bars)) {
      ht <- b$ht
      if (isTRUE(input$h1_blocknames)) ht <- tryCatch({ ht@ht_list[[1]]@column_title <- paste0("B", seq_len(b$nb)); ht }, error = function(e) b$ht)
      return(function() draw(ht))
    }
    mcs <- b$mc_order_ht; ro <- b$row_ord
    mat <- t(b$Z[mcs, ro, drop = FALSE]); rownames(mat) <- b$mods$name[ro]
    vals <- lapply(bars, function(k) mv(k, mcs))
    cols <- lapply(seq_along(bars), function(i) { v <- vals[[i]]
      if (is.numeric(v)) colorRamp2(quantile(v, c(.02, .5, .98), na.rm = TRUE) + c(0, 1e-9, 2e-9), c("#F3F6FB", "#7FA7D9", "#1F3F80")) else pal_for(level_order(bars[i], v)) })
    nms <- vapply(bars, key_label, ""); names(cols) <- nms
    adf <- as.data.frame(setNames(vals, nms), check.names = FALSE, stringsAsFactors = FALSE)
    top <- HeatmapAnnotation(df = adf, col = cols, annotation_name_gp = gpar(fontsize = 8), simple_anno_size = unit(3.2, "mm"),
                             annotation_legend_param = list(labels_gp = gpar(fontsize = 7), title_gp = gpar(fontsize = 8, fontface = "bold")))
    ht <- Heatmap(mat, name = "module z", col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick")),
                  cluster_rows = FALSE, cluster_columns = FALSE, column_split = factor(b$block[mcs], levels = seq_len(b$nb)),
                  cluster_column_slices = FALSE, column_title = if (isTRUE(input$h1_blocknames)) paste0("B", seq_len(b$nb)) else NULL,
                  column_title_gp = gpar(fontsize = 8), column_gap = unit(0.4, "mm"), show_column_names = FALSE,
                  row_names_gp = gpar(fontsize = min(8, 420 / nrow(mat))), top_annotation = top, use_raster = TRUE)
    function() draw(ht, merge_legend = TRUE)
  })
  output$h1_plot <- renderPlot({ f <- h1_draw(); f(); assign("01_module_heatmap", list(draw = f, w = input[["h1dl-w"]], h = input[["h1dl-h"]]), envir = figs) },
                               height = function() input$h1_height, res = 96)
  dl_server("h1dl", h1_draw, function() paste0(D()$name, "_module_heatmap"))

  # ============================================================ HEATMAP + GENES
  gene_list <- reactive({
    b <- req(D()); src <- input$g_src %||% "top"
    switch(src,
      top = unique(unlist(lapply(blocks_shown(), function(k) head(b$block_top[[k]], input$g_top_n %||% 2)))),
      var = { v <- b$bs$gene_var[b$bs$gene_frac >= 5e-5]; head(names(sort(v, decreasing = TRUE)), input$g_var_n %||% 40) },
      custom = match_genes(parse_list(input$g_custom), rownames(b$counts))$found,
      orig = intersect(ORIGINAL_LIST, rownames(b$counts)))
  })
  output$g_status <- renderUI({
    b <- req(D()); g <- gene_list()
    miss <- switch(input$g_src, custom = match_genes(parse_list(input$g_custom), rownames(b$counts))$missing,
                   orig = setdiff(ORIGINAL_LIST, rownames(b$counts)), character(0))
    hint(paste0(length(g), " genes"), if (length(miss)) paste0(" | not found: ", paste(head(miss, 12), collapse = ", "), if (length(miss) > 12) " ..."))
  })
  add_genes <- function(g) {
    cur <- if (input$g_src == "custom") parse_list(input$g_custom) else gene_list()
    updateRadioButtons(session, "g_src", selected = "custom")
    updateTextAreaInput(session, "g_custom", value = paste(unique(c(cur, g)), collapse = ", "))
    showNotification(paste0("Added ", length(g), " gene(s) to 'Heatmap + genes'. Click 'Update heatmap' there."), type = "message")
  }

  mc_order_file <- reactive({ req(input$mc_order_file); as.character(utils::read.csv(input$mc_order_file$datapath, header = FALSE, stringsAsFactors = FALSE)[, 1]) })
  observe({ disable("h2_custom_order", is.null(input$mc_order_file)) })

  h2_data <- eventReactive(list(input$h2_go, D()), {
    b <- req(D()); mcs <- mcs_shown()
    if (isTRUE(input$h2_custom_order) && !is.null(input$mc_order_file)) {
      mcs <- intersect(mc_order_file(), mcs)
      validate(need(length(mcs) > 0, "None of the names in the custom order file match the blocks shown."))
    }
    ro <- b$row_ord; M <- t(b$Z[mcs, ro, drop = FALSE]); rownames(M) <- b$mods$name[ro]
    g <- gene_list(); G <- NULL
    if (length(g)) {
      G <- as.matrix(b$counts[g, mcs, drop = FALSE]); G <- log2(1 + sweep(G, 2, b$libsize[mcs], "/") * 1e5)
      if (isTRUE(input$h2_genez)) { G <- t(scale(t(G))); G[is.na(G)] <- 0 }
    }
    bar <- input$h2_bar %||% "none"
    list(M = M, G = G, mcs = mcs, blk = b$block[mcs], bar = if (bar != "none") mv(bar, mcs) else NULL, bar_key = bar,
         genez = isTRUE(input$h2_genez), prop = input$h2_gprop, lab = input$h2_lab, sep = isTRUE(input$h2_sep))
  }, ignoreNULL = FALSE)

  h2_plot <- reactive({
    d <- h2_data(); n <- length(d$mcs)
    r <- rle(d$blk); ends <- cumsum(r$lengths); mids <- ends - r$lengths / 2 + 0.5
    show_sep <- d$sep && length(r$lengths) <= 300
    base <- theme_minimal(base_size = 11) + theme(axis.title = element_blank(), panel.grid = element_blank(),
      axis.text.y = element_text(size = d$lab, colour = "#3A3F4B"), axis.text.x = element_blank(),
      legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8), plot.margin = margin(2, 8, 2, 8))
    vl <- if (show_sep && length(ends) > 1) geom_vline(xintercept = ends[-length(ends)] + 0.5, linewidth = 0.25, colour = "#1D2129") else NULL
    long <- function(m) data.frame(x = rep(seq_len(ncol(m)), each = nrow(m)), y = factor(rep(rownames(m), ncol(m)), levels = rev(rownames(m))), v = as.vector(m))
    xs <- function(top = FALSE) if (top && show_sep) scale_x_continuous(expand = c(0, 0), breaks = mids, labels = paste0("B", r$values), position = "top") else scale_x_continuous(expand = c(0, 0))
    parts <- list(); hts <- numeric(0)
    first_top <- TRUE
    if (!is.null(d$bar)) {
      lv <- level_order(d$bar_key, d$bar)
      p0 <- ggplot(data.frame(x = seq_len(n), y = key_label(d$bar_key), v = factor(d$bar, levels = lv)), aes(x, y, fill = v)) + geom_raster() +
        xs(TRUE) + base + theme(axis.text.x.top = element_text(size = 8, colour = "#3A3F4B")) + scale_fill_manual(values = pal_for(lv), name = key_label(d$bar_key), na.value = GREY)
      parts <- c(parts, list(p0)); hts <- c(hts, 0.035); first_top <- FALSE
    }
    p1 <- ggplot(long(d$M), aes(x, y, fill = v)) + geom_raster() + xs(first_top) + vl + base +
      theme(axis.text.x.top = element_text(size = 8, colour = "#3A3F4B")) +
      scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", limits = c(-2, 2), oob = scales::squish, name = "Module\nscore (z)")
    gp <- if (!is.null(d$G) && d$prop > 0) d$prop else 0
    parts <- c(parts, list(p1)); hts <- c(hts, 1 - gp)
    if (gp > 0) {
      p2 <- ggplot(long(d$G), aes(x, y, fill = v)) + geom_raster() + xs() + vl + base +
        scale_fill_viridis_c(name = if (d$genez) "Gene\n(z)" else "Expression\nlog2(1+CPM/10)")
      parts <- c(parts, list(p2)); hts <- c(hts, gp)
    }
    wrap_plots(parts, ncol = 1, heights = hts) + plot_layout(guides = "collect")
  })
  output$h2_out <- renderUI({
    if (isTRUE(input$h2_hover)) plotly::plotlyOutput("h2_plotly", height = paste0(input$h2_height, "px"))
    else plotOutput("h2_plot", height = "auto")
  })
  output$h2_plot <- renderPlot({
    p <- h2_plot(); print(p)
    assign("02_heatmap_genes", list(draw = function() print(p), w = input[["h2dl-w"]], h = input[["h2dl-h"]]), envir = figs)
  }, height = function() input$h2_height, res = 96)
  output$h2_plotly <- plotly::renderPlotly({
    d <- h2_data(); xl <- paste0(d$mcs, " | B", d$blk)
    p1 <- plotly::plot_ly(x = xl, y = rownames(d$M), z = d$M, type = "heatmap", zmin = -2, zmax = 2,
                          colorscale = list(c(0, "navy"), c(0.5, "white"), c(1, "firebrick")), colorbar = list(title = "Module z", y = 0.8, len = 0.35))
    if (is.null(d$G) || d$prop == 0) return(plotly::layout(p1, xaxis = list(showticklabels = FALSE), yaxis = list(autorange = "reversed")))
    p2 <- plotly::plot_ly(x = xl, y = rownames(d$G), z = d$G, type = "heatmap", colorscale = "Viridis", colorbar = list(title = "Expression", y = 0.2, len = 0.35))
    plotly::layout(plotly::subplot(p1, p2, nrows = 2, heights = c(1 - d$prop, d$prop), shareX = TRUE, margin = 0.01),
                   xaxis = list(showticklabels = FALSE), yaxis = list(autorange = "reversed"), yaxis2 = list(autorange = "reversed"))
  })
  dl_server("h2dl", function() { p <- h2_plot(); function() print(p) }, function() paste0(D()$name, "_heatmap_genes"))

  # ============================================================ MODULE BY BLOCK
  output$bb_genes <- renderUI({ b <- req(D()); i <- match(req(input$bb_mod), b$mods$id); div(class = "mcx-genes", gsub(",", ", ", b$mods$genes[i])) })
  bb_plot <- reactive({
    b <- req(D()); id <- req(input$bb_mod); ks <- blocks_shown(); mcs <- mcs_shown()
    df <- data.frame(block = factor(paste0("B", b$block[mcs]), levels = paste0("B", ks)), z = b$Z[mcs, id])
    nm <- b$mods$name[match(id, b$mods$id)]
    p <- ggplot(df, aes(block, z))
    fill <- input$bb_fill %||% "none"
    if (fill != "none" && !is.null(mv(fill, mcs[1]))) {
      bm <- block_major(fill); df$f <- bm[b$block[mcs]]; lv <- level_order(fill, df$f)
      p <- ggplot(df, aes(block, z)) + geom_boxplot(aes(fill = factor(f, levels = lv)), width = 0.75, linewidth = 0.35, outlier.shape = NA, colour = "#3A3F4B") +
        scale_fill_manual(values = scales::alpha(pal_for(lv), 0.75), name = paste("Main", tolower(key_label(fill)), "of block"), na.value = GREY)
    } else p <- p + geom_boxplot(width = 0.75, fill = "#E8ECF8", colour = "#3A3F4B", linewidth = 0.35, outlier.shape = NA)
    p + geom_jitter(height = 0, width = 0.2, size = 0.35, alpha = 0.35, colour = "#1D2129") +
      geom_hline(yintercept = 0, linetype = 3, colour = "#9AA0AC") +
      labs(title = nm, x = NULL, y = "Module score (z)") + theme_mcx() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  })
  output$bb_plot <- renderPlot({ p <- bb_plot(); print(p); assign("03_module_by_block", list(draw = function() print(p), w = input[["bbdl-w"]], h = input[["bbdl-h"]]), envir = figs) },
                               height = function() input$bb_height, res = 96)
  dl_server("bbdl", function() { p <- bb_plot(); function() print(p) }, function() paste0(D()$name, "_", gsub(" ", "", D()$mods$name[match(input$bb_mod, D()$mods$id)]), "_by_block"))

  # ============================================================ ANNOTATE
  blk_summary <- reactive({
    b <- req(D()); cl <- CELLS(); a <- rv$ann; pl <- primary_label()
    annb <- vapply(seq_len(b$nb), function(k) { t <- sort(table(a$annotation[b$block == k]), decreasing = TRUE)
      if (length(t) == 1) names(t) else paste0(names(t)[1], " (", round(100 * t[1] / sum(t)), "%)") }, "")
    topm <- apply(b$block_z, 1, function(r) paste(b$mods$name[match(colnames(b$block_z)[order(-r)[1:3]], b$mods$id)], collapse = ", "))
    data.frame(Block = seq_len(b$nb), Metacells = b$bs$sizes,
               Cells = if (!is.null(cl$n_cells)) as.vector(tapply(cl$n_cells, b$block, sum)) else NA,
               `Cell type` = if (!is.null(pl)) block_major(paste0("lab:", pl)) else NA,
               `Your annotation` = annb, `Top modules` = topm,
               `Top genes` = vapply(b$block_top, function(g) paste(head(g, 4), collapse = ", "), ""),
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  output$sel_blk_tbl <- renderDT({
    s <- blk_summary(); s <- s[, colSums(!is.na(s)) > 0, drop = FALSE]
    datatable(s, rownames = FALSE, selection = list(mode = "multiple"),
              options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollY = "380px", scrollCollapse = TRUE, scrollX = TRUE))
  })
  blk_proxy <- dataTableProxy("sel_blk_tbl")
  observeEvent(input$sel_all, selectRows(blk_proxy, seq_len(D()$nb)))
  observeEvent(input$sel_none, selectRows(blk_proxy, NULL))
  observeEvent(input$sel_type_col, {
    req(D(), nzchar(input$sel_type_col)); updateSelectizeInput(session, "sel_type_val", choices = level_order(input$sel_type_col, mv(input$sel_type_col)), selected = character(0))
  })

  mc_range <- function(a, z) { b <- D(); i <- match(a, b$mc_order_ht); j <- match(z, b$mc_order_ht); if (is.na(i) || is.na(j)) character(0) else b$mc_order_ht[i:j] }
  sel_mcs <- reactive({
    b <- req(D())
    s <- switch(input$sel_mode %||% "block",
      block = b$mcs[b$block %in% input$sel_blk_tbl_rows_selected],
      type = if (length(input$sel_type_val)) b$mcs[mv(input$sel_type_col) %in% input$sel_type_val] else character(0),
      range = mc_range(input$sel_from, input$sel_to),
      names = intersect(parse_list(input$sel_names), b$mcs))
    setdiff(s, mc_range(input$ex_from, input$ex_to))
  })
  output$sel_info <- renderUI({
    b <- req(D()); s <- sel_mcs()
    if (!length(s)) return(div(class = "mcx-selinfo", icon("hand-pointer"), " Nothing selected yet."))
    ks <- sort(unique(b$block[s]))
    extra <- if (input$sel_mode == "names") { bad <- setdiff(parse_list(input$sel_names), b$mcs); if (length(bad)) paste0(" | not found: ", paste(head(bad, 8), collapse = ", ")) }
    div(class = "mcx-selinfo", icon("circle-check"), strong(paste0(" ", fmt_n(length(s)), " metacells")),
        paste0(" from ", if (length(ks) > 12) paste(length(ks), "blocks") else paste0("block", if (length(ks) > 1) "s", " ", paste(ks, collapse = ", "))), extra)
  })

  autosave_file <- reactive(paste0(AUTOSAVE_DIR, "/", gsub("[^A-Za-z0-9_]+", "_", req(D())$name), "_annotation.csv"))
  ann_export <- function() {
    b <- D(); cl <- CELLS(); a <- rv$ann
    out <- data.frame(metacell = a$metacell, annotation = a$annotation, block = b$block[a$metacell], stringsAsFactors = FALSE)
    if (!is.null(cl$n_cells)) out$n_cells <- cl$n_cells[match(a$metacell, b$mcs)]
    for (c in names(cl$labels)) out[[c]] <- cl$labels[[c]][match(a$metacell, b$mcs)]
    out
  }
  log_it <- function(action, n, name) rv$log <- rbind(data.frame(time = format(Sys.time(), "%H:%M:%S"), action = action, metacells = as.integer(n), name = name), rv$log)
  set_ann <- function(new, action, n, name) {
    rv$hist <- c(list(rv$ann), head(rv$hist, 29)); rv$ann <- new; log_it(action, n, name)
    if (isTRUE(input$autosave)) tryCatch({ dir.create(AUTOSAVE_DIR, recursive = TRUE, showWarnings = FALSE); utils::write.csv(ann_export(), autosave_file(), row.names = FALSE) },
                                        error = function(e) showNotification("Autosave failed (folder not writable).", type = "warning"))
  }
  observeEvent(input$ann_assign, {
    s <- sel_mcs(); nm <- trimws(input$ann_name)
    if (!length(s)) return(showNotification("Select some metacells first.", type = "warning"))
    if (!nzchar(nm)) return(showNotification("Type a name first.", type = "warning"))
    a <- rv$ann; a$annotation[a$metacell %in% s] <- nm; set_ann(a, "named", length(s), nm)
    selectRows(blk_proxy, NULL); updateTextInput(session, "ann_name", value = "")
    showNotification(HTML(paste0("<b>", fmt_n(length(s)), "</b> metacells named <b>", htmltools::htmlEscape(nm), "</b>")), type = "message")
  })
  observeEvent(input$ann_undo, {
    if (!length(rv$hist)) return(showNotification("Nothing to undo.", type = "warning"))
    rv$ann <- rv$hist[[1]]; rv$hist <- rv$hist[-1]; log_it("undo", NA, "")
  })
  observeEvent(input$ann_init, {
    req(input$ann_init_col != "none"); b <- D(); a <- rv$ann
    v <- mv(input$ann_init_col, a$metacell); v[is.na(v)] <- "undefined"; a$annotation <- v
    set_ann(a, paste("started from", key_label(input$ann_init_col)), nrow(a), "")
  })
  load_ann <- function(path, label) {
    x <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(x) || !all(c("metacell", "annotation") %in% names(x))) { rv$ann_msg <- "The file needs the columns 'metacell' and 'annotation'."; return() }
    a <- rv$ann; i <- match(a$metacell, x$metacell); a$annotation[!is.na(i)] <- as.character(x$annotation[i[!is.na(i)]])
    set_ann(a, paste("loaded", label), sum(!is.na(i)), "")
    rv$ann_msg <- paste0(fmt_n(sum(!is.na(i))), " metacells matched", if (any(!x$metacell %in% a$metacell)) paste0("; ", sum(!x$metacell %in% a$metacell), " names in the file are not in this dataset") else "", ".")
  }
  observeEvent(input$ann_file, load_ann(input$ann_file$datapath, input$ann_file$name))
  output$ann_load_msg <- renderUI({ m <- rv$ann_msg; if (is.null(m) || m == "autosave") NULL else hint(m) })
  output$autosave_ui <- renderUI({
    f <- autosave_file()
    tagList(hint("Saved to ", code(f)),
            if (identical(rv$ann_msg, "autosave") && file.exists(f)) tagList(
              div(class = "alert alert-info py-2 px-3 small", "An autosaved annotation for this dataset exists (", format(file.mtime(f), "%d %b %H:%M"), ")."),
              actionButton("ann_restore", "Restore it", class = "btn-sm btn-outline-primary")))
  })
  observeEvent(input$ann_restore, { load_ann(autosave_file(), "autosave"); rv$ann_msg <- NULL })

  output$ann_tbl <- renderDT({
    b <- req(D()); a <- req(rv$ann); cl <- CELLS()
    t <- as.data.frame(table(a$annotation), stringsAsFactors = FALSE); names(t) <- c("Name", "Metacells")
    t <- t[order(t$Name %in% NO_LABEL, -t$Metacells), ]
    if (!is.null(cl$n_cells)) t$Cells <- as.integer(tapply(cl$n_cells[match(a$metacell, b$mcs)], a$annotation, sum)[t$Name])
    t$`%` <- round(100 * t$Metacells / sum(t$Metacells), 1)
    t$Name <- paste0("<span class='swatch' style='background:", pal_for(t$Name), "'></span>", htmltools::htmlEscape(t$Name))
    datatable(t, rownames = FALSE, escape = FALSE, selection = "none", options = list(dom = "t", paging = FALSE, ordering = FALSE))
  })
  output$ann_log <- renderDT(datatable(rv$log %||% data.frame(), rownames = FALSE, selection = "none", options = list(dom = "tp", pageLength = 5, ordering = FALSE)))

  relabel_code <- function() {
    b <- D(); a <- ann_export()
    bl <- do.call(rbind, lapply(split(a, a$block), function(d) { t <- sort(table(d$annotation), decreasing = TRUE)
      data.frame(block = d$block[1], label = names(t)[1], purity = as.numeric(t[1]) / nrow(d)) }))
    bl <- bl[order(bl$block), ]; named <- bl[!bl$label %in% NO_LABEL, ]
    dev <- a[!a$annotation %in% NO_LABEL & a$annotation != bl$label[match(a$block, bl$block)], ]
    q <- function(x) paste0("\"", x, "\"")
    c(paste0("# MCexplorer | ", b$name, " | ", format(Sys.time(), "%Y-%m-%d %H:%M")),
      "# Block numbers are left-to-right heatmap positions (same as 'block' in the modules scoring script).",
      "# Each block gets the name most of its metacells received:",
      if (nrow(named)) paste0("BLOCK_RELABEL <- c(", paste0("\n  ", q(named$block), " = ", q(named$label), collapse = ","), "\n)") else "BLOCK_RELABEL <- c()",
      if (nrow(named)) paste0("# share of the block with that name: ", paste0(named$block, " = ", round(named$purity, 2), collapse = ", ")),
      "", "# Metacells named differently from the rest of their block:",
      if (nrow(dev)) paste0("MC_RELABEL <- c(", paste0("\n  ", q(dev$metacell), " = ", q(dev$annotation), collapse = ","), "\n)") else "MC_RELABEL <- c()")
  }
  block_ann_table <- function() { a <- ann_export(); t <- as.data.frame.matrix(table(a$block, a$annotation)); cbind(block = as.integer(rownames(t)), t) }
  output$ann_dl_csv <- downloadHandler(function() paste0(D()$name, "_annotation_", Sys.Date(), ".csv"), function(file) utils::write.csv(ann_export(), file, row.names = FALSE))
  output$ann_dl_r   <- downloadHandler(function() paste0(D()$name, "_relabel_", Sys.Date(), ".R"), function(file) writeLines(relabel_code(), file))
  output$ann_dl_blk <- downloadHandler(function() paste0(D()$name, "_blocks_x_annotation_", Sys.Date(), ".csv"), function(file) utils::write.csv(block_ann_table(), file, row.names = FALSE))

  # ============================================================ COMPARE
  cg_levels <- reactive({ req(D(), input$cg_by); v <- mv(input$cg_by); if (is.null(v)) character(0) else level_order(input$cg_by, v) })
  observeEvent(cg_levels(), {
    lv <- cg_levels(); cur <- isolate(input$cg_groups)
    keep <- intersect(cur, lv); sel <- if (length(keep) >= 2) c(keep, setdiff(lv, isolate(rv$cg_prev) %||% character(0))) else lv
    rv$cg_prev <- lv
    updateSelectizeInput(session, "cg_groups", choices = lv, selected = intersect(sel, lv))
  })
  observeEvent(input$cg_groups, {
    g <- input$cg_groups
    if (length(g) > 4 && isolate(input$cg_cmp) == "all") updateRadioButtons(session, "cg_cmp", selected = "ref")
    updateSelectInput(session, "cg_ref", choices = g, selected = if (isTRUE(isolate(input$cg_ref) %in% g)) isolate(input$cg_ref) else g[1])
    pr <- if (length(g) > 1) vapply(combn(g, 2, simplify = FALSE), paste, "", collapse = " vs ") else character(0)
    updateSelectizeInput(session, "cg_pairs", choices = pr, selected = intersect(isolate(input$cg_pairs), pr))
  }, ignoreNULL = FALSE)
  observe({
    ng <- length(input$cg_groups); st <- input$cg_stats %||% "none"
    disable("cg_stats", ng < 2)
    off <- ng < 2 || st == "none"
    disable(c("cg_cmp", "cg_adj", "cg_hidens"), off)
    disable("cg_ref", off || input$cg_cmp != "ref")
    disable("cg_pairs", off || input$cg_cmp != "pick")
  })
  observeEvent(input$cg_stats, update_switch("cg_hidens", value = input$cg_stats == "stars"), ignoreInit = TRUE)

  cg_result <- reactive({
    b <- req(D()); key <- req(input$cg_by); groups <- input$cg_groups
    validate(need(length(groups) > 0, "Choose at least one group to show."))
    grp <- mv(key); validate(need(!is.null(grp), "These groups need cell information.")); sel <- b$mcs[!is.na(grp) & grp %in% groups]
    what <- input$cg_what; ylab <- ""; title <- ""
    if (what == "module") { id <- req(input$cg_mod); score <- b$Z[sel, id]; title <- b$mods$name[match(id, b$mods$id)]; ylab <- "Module score (z)" }
    else if (what == "gene") { g <- req(input$cg_gene); score <- log2(1 + as.vector(b$counts[g, sel]) / b$libsize[sel] * 1e5); title <- g; ylab <- "Expression, log2(1 + CPM/10)" }
    else if (what == "set") {
      m <- match_genes(parse_list(input$cg_set), rownames(b$counts))
      validate(need(length(m$found) >= 2, paste0("Type at least two genes.", if (length(m$missing)) paste0(" Not found: ", paste(m$missing, collapse = ", ")) else "")))
      v <- Matrix::colSums(log1p(b$counts[m$found, sel, drop = FALSE])) / b$allsum[sel]; score <- as.vector(scale(v))
      title <- paste0("Gene set (", length(m$found), " genes)"); ylab <- "Gene set score (z)"
    } else { k <- req(input$cg_metric); score <- mv(k, sel); validate(need(!is.null(score), "This metric needs cell information.")); title <- key_label(k); ylab <- key_label(k) }
    df <- data.frame(metacell = sel, group = factor(grp[match(sel, b$mcs)], levels = groups), score = score)
    df <- df[!is.na(df$score), ]
    pal <- pal_for(groups)
    p <- ggplot(df, aes(group, score)) +
      geom_boxplot(aes(fill = group, colour = group), width = 0.7, linewidth = 0.5, outlier.shape = NA, alpha = 0.35) +
      geom_jitter(aes(colour = group), height = 0, width = 0.18, size = input$cg_pt, alpha = 0.55) +
      scale_fill_manual(values = pal, drop = FALSE) + scale_colour_manual(values = pal, drop = FALSE) +
      labs(title = title, x = NULL, y = ylab, fill = key_label(key), colour = key_label(key)) + theme_mcx(input$cg_txt) +
      theme(axis.text.x = element_text(angle = if (max(nchar(groups)) > 6) 40 else 0, hjust = if (max(nchar(groups)) > 6) 1 else 0.5),
            legend.position = if (isTRUE(input$cg_legend)) "right" else "none")
    if (isTRUE(input$cg_mean)) p <- p + geom_hline(yintercept = mean(df$score), linetype = 2, colour = "#9AA0AC", linewidth = 0.4)

    stats <- NULL
    if (input$cg_stats != "none" && length(groups) > 1) {
      ok <- names(which(table(df$group) >= 2)); dfs <- droplevels(df[df$group %in% ok, ])
      kw <- tryCatch(rstatix::kruskal_test(dfs, score ~ group), error = function(e) NULL)
      comps <- switch(input$cg_cmp,
        all = if (length(ok) > 1) combn(ok, 2, simplify = FALSE) else list(),
        ref = lapply(setdiff(ok, input$cg_ref), function(g) c(input$cg_ref, g)),
        pick = lapply(strsplit(input$cg_pairs %||% character(0), " vs ", fixed = TRUE), function(x) x))
      comps <- Filter(function(x) length(x) == 2 && all(x %in% ok), comps)
      adj_name <- c(bonferroni = "Bonferroni", BH = "Benjamini-Hochberg", holm = "Holm", none = "no")[[input$cg_adj]]
      if (!is.null(kw)) p <- p + labs(subtitle = paste0("Kruskal-Wallis p = ", signif(kw$p, 3), "   |   pairwise Wilcoxon, ", adj_name, " correction"))
      stats <- if (!is.null(kw)) data.frame(group1 = "all shown", group2 = "(Kruskal-Wallis)", n1 = nrow(dfs), n2 = NA, statistic = unname(kw$statistic), p = kw$p, p.adj = NA, signif = "")
      if (length(comps)) {
        pw <- tryCatch(rstatix::wilcox_test(dfs, score ~ group, comparisons = comps, p.adjust.method = input$cg_adj), error = function(e) NULL)
        if (!is.null(pw)) {
          if (!"p.adj" %in% names(pw)) pw$p.adj <- pw$p
          pw <- rstatix::add_significance(pw, "p.adj", output.col = "signif")
          pw$lab <- ifelse(pw$p.adj < 1e-3, formatC(pw$p.adj, format = "e", digits = 1), signif(pw$p.adj, 2))
          pw <- rstatix::add_xy_position(pw, x = "group", step.increase = 0.12)
          p <- p + ggpubr::stat_pvalue_manual(pw, label = if (input$cg_stats == "p") "lab" else "signif", tip.length = 0.01, vjust = -0.25,
                                              hide.ns = isTRUE(input$cg_hidens), size = input$cg_statsize, bracket.size = 0.35, inherit.aes = FALSE) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))
          pwd <- data.frame(group1 = pw$group1, group2 = pw$group2, n1 = pw$n1, n2 = pw$n2, statistic = pw$statistic, p = pw$p, p.adj = pw$p.adj, signif = pw$signif)
          stats <- rbind(stats, pwd)
        }
      }
    }
    list(plot = p, data = df, stats = stats, name = paste0(D()$name, "_", gsub("[^A-Za-z0-9]+", "_", title), "_by_", gsub("[^A-Za-z0-9]+", "_", key_label(key))))
  })
  output$cg_plot <- renderPlot({
    r <- cg_result(); print(r$plot)
    assign("05_compare", list(draw = function() print(r$plot), w = input[["cgdl-w"]], h = input[["cgdl-h"]]), envir = figs)
    if (!is.null(r$stats)) assign("compare_statistics", r$stats, envir = tables)
  }, height = function() input$cg_height, res = 96)
  dl_server("cgdl", function() { p <- cg_result()$plot; function() print(p) }, function() cg_result()$name)
  output$cg_stats_tbl <- renderDT({
    s <- cg_result()$stats
    if (is.null(s)) return(datatable(data.frame(` ` = "Choose 'Show p-values' or 'Show stars' to run statistics.", check.names = FALSE), rownames = FALSE, selection = "none", options = list(dom = "t")))
    datatable(s, rownames = FALSE, selection = "none", options = list(dom = "tp", pageLength = 8)) |> formatSignif(c("statistic", "p", "p.adj"), 3)
  })
  output$cg_dl_stats <- downloadHandler(function() paste0(cg_result()$name, "_stats.csv"), function(file) utils::write.csv(cg_result()$stats %||% data.frame(), file, row.names = FALSE))

  # colour pickers (shared colour memory, so a name keeps its colour in every plot)
  picker_ui <- function(prefix, lv) {
    if (!length(lv)) return(hint("No groups."))
    if (length(lv) > 30) return(hint("Too many groups to colour by hand (", length(lv), ")."))
    pal <- isolate(pal_for(lv))
    tagList(lapply(seq_along(lv), function(i) div(class = "mcx-colrow",
      tags$input(type = "color", id = paste0(prefix, i), class = "mcx-color", value = substr(grDevices::rgb(t(grDevices::col2rgb(pal[[i]])), maxColorValue = 255), 1, 7)),
      span(lv[i]))))
  }
  picker_watch <- function(prefix, lv_reactive) observe({
    lv <- lv_reactive(); vals <- lapply(seq_along(lv), function(i) input[[paste0(prefix, i)]])
    if (!length(lv) || any(vapply(vals, is.null, TRUE))) return()
    vals <- unlist(vals); cur <- isolate(pal_for(lv))
    ch <- toupper(vals) != toupper(substr(cur, 1, 7))
    if (any(ch)) { new <- isolate(rv$colors); new[lv[ch]] <- vals[ch]; rv$colors <- new }
  })
  output$cg_colors <- renderUI({ rv$nonce; picker_ui("cgcol_", input$cg_groups) })
  picker_watch("cgcol_", reactive(input$cg_groups))
  observeEvent(input$cg_colreset, { g <- input$cg_groups; rv$colors <- rv$colors[setdiff(names(rv$colors), g)]; rv$nonce <- rv$nonce + 1 })

  # ============================================================ COMPOSITION
  cp_levels <- reactive({ cl <- req(CELLS()); col <- req(input$cp_col); req(col %in% names(cl$comp)); colnames(cl$comp[[col]]) })
  output$cp_colors <- renderUI(picker_ui("cpcol_", cp_levels()))
  picker_watch("cpcol_", cp_levels)
  cp_result <- reactive({
    b <- req(D()); cl <- req(CELLS()); col <- req(input$cp_col); by <- req(input$cp_by); req(col %in% names(cl$comp))
    cm <- cl$comp[[col]]
    g <- mv(by); lv_g <- level_order(by, g); keep <- !is.na(g) & g %in% lv_g
    agg <- rowsum(cm[keep, , drop = FALSE], g[keep]); ncell <- rowSums(agg)
    if (input$cp_scale == "norm") agg <- sweep(agg, 2, pmax(colSums(cm), 1), "/")
    frac <- agg / pmax(rowSums(agg), 1e-12)
    df <- data.frame(group = factor(rep(rownames(frac), ncol(frac)), levels = intersect(lv_g, rownames(frac))),
                     level = factor(rep(colnames(frac), each = nrow(frac)), levels = colnames(cm)), share = as.vector(frac),
                     cells_in_group = as.vector(ncell[rownames(frac)]))
    p <- ggplot(df, aes(group, share, fill = level)) + geom_col(width = 0.82, colour = "white", linewidth = 0.15) +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) + scale_fill_manual(values = pal_for(colnames(cm)), name = pretty_name(col)) +
      labs(title = paste(pretty_name(col), "by", tolower(key_label(by))), x = NULL,
           y = if (input$cp_scale == "norm") "Share (adjusted for group size)" else "Share of cells") + theme_mcx() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.y = element_blank())
    if (isTRUE(input$cp_flip)) p <- p + coord_flip()
    list(plot = p, data = df)
  })
  output$cp_plot <- renderPlot({
    r <- cp_result(); print(r$plot)
    assign("06_composition", list(draw = function() print(r$plot), w = input[["cpdl-w"]], h = input[["cpdl-h"]]), envir = figs)
    assign("composition", r$data, envir = tables)
  }, height = function() input$cp_height, res = 96)
  dl_server("cpdl", function() { p <- cp_result()$plot; function() print(p) }, function() paste0(D()$name, "_composition_", input$cp_col))

  # ============================================================ FIND GENES
  output$fg_info <- renderUI({
    b <- req(D()); g <- req(input$fg_gene); m <- b$gene2mod[g]
    if (is.na(m)) return(div(class = "mcx-genes", strong(g), " is not part of any gene module (it was not variable enough to be used for modules)."))
    i <- match(m, b$mods$name)
    tagList(div(class = "mcx-section-title", paste0(m, " (", b$mods$n_genes[i], " genes)")), div(class = "mcx-genes", gsub(",", ", ", b$mods$genes[i])))
  })
  fg_plot <- reactive({
    b <- req(D()); g <- req(input$fg_gene); mcs <- mcs_shown(); ks <- blocks_shown(); pl <- primary_label()
    df <- data.frame(block = factor(paste0("B", b$block[mcs]), levels = paste0("B", ks)), v = log2(1 + as.vector(b$counts[g, mcs]) / b$libsize[mcs] * 1e5))
    if (!is.null(pl)) { df$f <- block_major(paste0("lab:", pl))[b$block[mcs]]; lv <- level_order(paste0("lab:", pl), df$f) }
    p <- ggplot(df, aes(block, v)) +
      (if (!is.null(pl)) geom_boxplot(aes(fill = factor(f, levels = lv)), outlier.shape = NA, width = 0.75, linewidth = 0.3, colour = "#3A3F4B")
       else geom_boxplot(fill = "#E8ECF8", outlier.shape = NA, width = 0.75, linewidth = 0.3, colour = "#3A3F4B")) +
      geom_jitter(height = 0, width = 0.2, size = 0.3, alpha = 0.3) +
      labs(title = g, x = NULL, y = "Expression, log2(1 + CPM/10)") + theme_mcx() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    if (!is.null(pl)) p <- p + scale_fill_manual(values = scales::alpha(pal_for(lv), 0.75), name = paste("Main", pretty_name(pl), "of block"), na.value = GREY)
    p
  })
  output$fg_plot <- renderPlot({ p <- fg_plot(); print(p); assign("07_gene_by_block", list(draw = function() print(p), w = input[["fgdl-w"]], h = input[["fgdl-h"]]), envir = figs) }, height = 420, res = 96)
  dl_server("fgdl", function() { p <- fg_plot(); function() print(p) }, function() paste0(D()$name, "_", input$fg_gene, "_by_block"))
  observeEvent(input$fg_add, add_genes(input$fg_gene))

  observe(disable("tg_b", input$tg_vs != "blocks"))
  tg_result <- reactive({
    b <- req(D()); A <- as.integer(sub("^B", "", req(input$tg_a)))
    B <- if (input$tg_vs == "blocks") setdiff(as.integer(sub("^B", "", input$tg_b %||% character(0))), A) else NULL
    validate(need(input$tg_vs == "rest" || length(B), "Choose blocks to compare with."))
    r <- head(rank_genes(b$bs, A, B), input$tg_n %||% 50)
    r$module <- unname(b$gene2mod[r$gene]); r$log2FC <- round(r$log2FC, 2); r$mean_frac <- signif(r$mean_frac, 3)
    assign("top_genes", r, envir = tables); r
  })
  output$tg_tbl <- renderDT(datatable({ validate(need(length(input$tg_a) > 0, "Pick one or more blocks.")); tg_result() }, rownames = FALSE, selection = "none", colnames = c("Gene", "log2 fold change", "Mean fraction", "Module"), options = list(pageLength = 15)))
  output$tg_dl <- downloadHandler(function() paste0(D()$name, "_top_genes_", paste(input$tg_a, collapse = "_"), ".csv"), function(file) utils::write.csv(tg_result(), file, row.names = FALSE))
  observeEvent(input$tg_add, add_genes(head(tg_result()$gene, 20)))

  output$tm_tbl <- renderDT({
    b <- req(D()); k <- as.integer(sub("^B", "", req(input$tm_block))); z <- b$block_z[k, ]
    o <- order(z, decreasing = TRUE); o <- o[z[o] > input$tm_thr]
    i <- match(names(z)[o], b$mods$id)
    datatable(data.frame(Module = b$mods$name[i], `Average score` = round(z[o], 2), Genes = vapply(b$mods$genes[i], function(g) paste(head(strsplit(g, ",")[[1]], 15), collapse = ", "), ""), check.names = FALSE),
              rownames = FALSE, selection = "none", options = list(pageLength = 15, dom = "tp"))
  })

  output$vg_note <- renderUI({ if (is.null(CELLS()$total_umis)) hint(icon("circle-info"), " Needs cell information with UMI counts (nCount_RNA or total_umis). Add it on the Data page.") })
  observeEvent(input$vg_go, {
    b <- req(D()); cl <- CELLS(); ks <- as.integer(sub("^B", "", input$vg_blocks %||% character(0)))
    if (!length(ks) || is.null(cl$total_umis)) return()
    mcs <- b$mcs[b$block %in% ks]; if (length(mcs) < 3) return(showNotification("Need at least 3 metacells.", type = "warning"))
    withProgress(message = "Testing genes", {
      umi <- cl$total_umis[match(mcs, b$mcs)]
      O <- as(b$counts[, mcs, drop = FALSE] %*% Diagonal(x = umi / b$libsize[mcs]), "CsparseMatrix")
      rs <- Matrix::rowSums(O); mx <- numeric(nrow(O))
      if (length(O@x)) { v <- tapply(O@x, O@i + 1L, max); mx[as.integer(names(v))] <- v }
      keep <- mx > 3 & rs >= 10; O <- O[keep, , drop = FALSE]; rs <- rs[keep]
      p <- umi / sum(umi); O2 <- O; O2@x <- O2@x^2
      chi <- as.vector(O2 %*% (1 / p)) / rs - rs
      pv <- stats::pchisq(chi, df = length(mcs) - 1, lower.tail = FALSE)
      res <- data.frame(gene = rownames(O), chi2 = round(chi, 1), p = signif(pv, 3), p_adj = signif(stats::p.adjust(pv, "BH"), 3),
                        umis = round(rs), module = unname(b$gene2mod[rownames(O)]), stringsAsFactors = FALSE)
      res <- res[order(-res$chi2), ]
      rv$vg <- res; assign("genes_varying_within_blocks", res, envir = tables)
    })
  })
  output$vg_tbl <- renderDT({ req(rv$vg); datatable(head(rv$vg, 1000), rownames = FALSE, selection = "none", colnames = c("Gene", "Chi-square", "p", "Adjusted p (BH)", "UMIs", "Module"), options = list(pageLength = 15)) })
  output$vg_dl <- downloadHandler(function() paste0(D()$name, "_genes_varying_within_blocks.csv"), function(file) utils::write.csv(req(rv$vg), file, row.names = FALSE))
  observeEvent(input$vg_add, { req(rv$vg); add_genes(head(rv$vg$gene[rv$vg$p_adj < 0.01], 30)) })

  mod_table <- reactive({
    b <- req(D()); o <- order(b$mods$num)
    top3 <- apply(b$block_z, 2, function(z) { k <- order(z, decreasing = TRUE)[1:min(3, length(z))]; paste0("B", k, " (", round(z[k], 1), ")", collapse = ", ") })
    data.frame(Module = b$mods$name[o], Genes = b$mods$n_genes[o], `Highest in blocks` = top3[b$mods$id[o]], `Gene list` = gsub(",", ", ", b$mods$genes[o]), check.names = FALSE)
  })
  output$mod_tbl <- renderDT(datatable(mod_table(), rownames = FALSE, selection = "none", options = list(pageLength = 15, scrollX = TRUE)))
  output$mod_dl <- downloadHandler(function() paste0(D()$name, "_gene_modules.csv"), function(file) utils::write.csv(mod_table(), file, row.names = FALSE))
  output$blk_tbl <- renderDT({ s <- blk_summary(); datatable(s[, colSums(!is.na(s)) > 0, drop = FALSE], rownames = FALSE, selection = "none", options = list(pageLength = 50, scrollX = TRUE)) })
  output$blk_dl <- downloadHandler(function() paste0(D()$name, "_blocks.csv"), function(file) utils::write.csv(blk_summary(), file, row.names = FALSE))

  # ============================================================ UMAP
  compute_umap <- function() {
    b <- req(D())
    if (!requireNamespace("uwot", quietly = TRUE)) { rv$umap_msg <- "The UMAP needs the 'uwot' package. Install it once with: install.packages(\"uwot\")"; return() }
    withProgress(message = "Computing UMAP from module scores", value = 0.3, {
      mcx_log("Computing UMAP (", nrow(b$Z), " metacells x ", ncol(b$Z), " module scores)")
      set.seed(42)
      xy <- uwot::umap(b$Z, n_neighbors = min(input$um_nn, nrow(b$Z) - 1), min_dist = input$um_md, n_threads = max(1, min(4, parallel::detectCores() - 1)))
      u <- list(xy = data.frame(x = xy[, 1], y = xy[, 2], row.names = b$mcs),
                source = sprintf("Computed from module scores (neighbours %d, minimum distance %.2f)", input$um_nn, input$um_md))
      UMAP(u); cache_write(b$cache_file, b$cache_key, "umap", u); rv$umap_msg <- NULL
    })
  }
  observeEvent(input$nav, { if (input$nav == "umap" && !is.null(D()) && is.null(UMAP())) compute_umap() })
  observeEvent(input$um_recompute, compute_umap())
  output$um_source <- renderUI(hint(UMAP()$source %||% rv$umap_msg %||% ""))
  um_plot <- reactive({
    b <- req(D()); u <- UMAP(); validate(need(!is.null(u), rv$umap_msg %||% "Computing UMAP..."))
    xy <- u$xy; mcs <- rownames(xy); by <- req(input$um_by); num <- FALSE
    if (by == "module") { v <- b$Z[mcs, req(input$um_mod)]; title <- b$mods$name[match(input$um_mod, b$mods$id)]; num <- TRUE }
    else if (by == "gene") { g <- req(input$um_gene); v <- log2(1 + as.vector(b$counts[g, mcs]) / b$libsize[mcs] * 1e5); title <- g; num <- TRUE }
    else { v <- mv(by, mcs); validate(need(!is.null(v), "This needs cell information.")); title <- key_label(by); num <- is.numeric(v) }
    xy$v <- v; if (num) xy <- xy[order(xy$v), ] else { lv <- level_order(by, v); xy$v <- factor(xy$v, levels = lv) }
    p <- ggplot(xy, aes(x, y, colour = v)) + geom_point(size = input$um_pt, stroke = 0) + coord_equal() +
      labs(title = title, x = "UMAP 1", y = "UMAP 2", colour = NULL) + theme_mcx() +
      theme(axis.text = element_blank(), panel.grid = element_blank(), axis.line.x = element_blank(),
            panel.border = element_rect(fill = NA, colour = "#E5E7EB"))
    if (num) p <- p + if (by == "module") scale_colour_gradient2(low = "navy", mid = "#E9EAEE", high = "firebrick", limits = c(-2, 2), oob = scales::squish) else scale_colour_viridis_c()
    else {
      p <- p + scale_colour_manual(values = pal_for(lv), na.value = GREY) + guides(colour = guide_legend(override.aes = list(size = 3.5)))
      if (isTRUE(input$um_labels) && length(lv) <= 60) {
        cen <- stats::aggregate(cbind(x, y) ~ v, data = xy, FUN = stats::median)
        p <- p + geom_label(data = cen, aes(x, y, label = v), inherit.aes = FALSE, size = 3, label.size = 0, fill = scales::alpha("white", 0.8), colour = "#1D2129")
      }
    }
    p
  })
  output$um_plot <- renderPlot({ p <- um_plot(); print(p); assign("08_umap", list(draw = function() print(p), w = input[["umdl-w"]], h = input[["umdl-h"]]), envir = figs) },
                               height = function() input$um_height, res = 96)
  dl_server("umdl", function() { p <- um_plot(); function() print(p) }, function() paste0(D()$name, "_umap_", gsub("[^A-Za-z0-9]+", "_", input$um_by)))

  # ============================================================ SAVE
  session_obj <- function() list(app = "MCexplorer", version = APP_VERSION, saved = Sys.time(), dataset = D()$name,
    ann = rv$ann, log = rv$log, colors = rv$colors, inputs = lapply(setNames(names(RESTORE), names(RESTORE)), function(id) isolate(input[[id]])))
  apply_session <- function(s) {
    for (id in names(s$inputs)) {
      v <- s$inputs[[id]]; type <- RESTORE[[id]]
      if (is.null(type) || (is.null(v) && type != "multi")) next
      switch(type, text = updateTextAreaInput(session, id, value = v), numeric = updateNumericInput(session, id, value = v),
             slider = updateSliderInput(session, id, value = v), switch = update_switch(id, value = v),
             radio = updateRadioButtons(session, id, selected = v), select = updateSelectInput(session, id, selected = v),
             multi = updateSelectizeInput(session, id, selected = v %||% character(0)))
    }
    a <- rv$ann; i <- match(a$metacell, s$ann$metacell); a$annotation[!is.na(i)] <- s$ann$annotation[i[!is.na(i)]]
    rv$ann <- a; rv$log <- rbind(data.frame(time = format(Sys.time(), "%H:%M:%S"), action = "session restored", metacells = sum(!is.na(i)), name = ""), s$log %||% NULL)
    rv$colors <- s$colors %||% character(0)
  }
  observeEvent(input$session_file, {
    s <- tryCatch(readRDS(input$session_file$datapath), error = function(e) NULL)
    if (is.null(s) || !identical(s$app, "MCexplorer")) { rv$session_msg <- "This is not an MCexplorer session file."; return() }
    if (!identical(s$dataset, D()$name)) { rv$session_msg <- paste0("This session belongs to '", s$dataset, "'. Open that dataset first."); return() }
    apply_session(s); rv$session_msg <- paste0("Restored the session saved on ", format(s$saved, "%d %b %Y %H:%M"), ".")
  })
  output$session_msg <- renderUI(if (!is.null(rv$session_msg)) div(class = "alert alert-info py-2 px-3 small mt-2", rv$session_msg))
  output$dl_session <- downloadHandler(function() paste0(D()$name, "_session_", format(Sys.time(), "%Y%m%d_%H%M"), ".rds"), function(file) saveRDS(session_obj(), file))

  output$dl_all <- downloadHandler(
    function() paste0(D()$name, "_MCexplorer_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip"),
    function(file) {
      b <- D(); td <- tempfile("mcx_"); for (s in c("annotation", "tables", "figures")) dir.create(paste0(td, "/", s), recursive = TRUE)
      wcsv <- function(x, f) utils::write.csv(x, paste0(td, "/", f), row.names = FALSE)
      withProgress(message = "Preparing zip", {
        wcsv(ann_export(), "annotation/metacell_annotation.csv"); wcsv(rv$log %||% data.frame(), "annotation/history.csv")
        wcsv(block_ann_table(), "annotation/blocks_x_annotation.csv"); writeLines(relabel_code(), paste0(td, "/annotation/relabel_code.R"))
        wcsv(mod_table(), "tables/gene_modules.csv"); wcsv(blk_summary(), "tables/blocks.csv")
        for (n in ls(tables)) wcsv(get(n, envir = tables), paste0("tables/", n, ".csv"))
        incProgress(0.3, detail = "figures")
        for (n in ls(figs)) { f <- get(n, envir = figs)
          tryCatch(save_plot(paste0(td, "/figures/", n, ".pdf"), f$draw, f$w, f$h, "pdf"), error = function(e) writeLines(conditionMessage(e), paste0(td, "/figures/", n, "_error.txt"))) }
        saveRDS(session_obj(), paste0(td, "/session.rds"))
        writeLines(c(paste("MCexplorer", APP_VERSION, "|", b$name), paste("Exported", format(Sys.time())), "",
                     "annotation/  your names per metacell, per block, the relabel code and the history",
                     "tables/      gene modules, blocks, and the last result of each table in the app",
                     "figures/     every figure looked at in this session, at the sizes set under each figure",
                     "session.rds  load it on the Save page to restore annotation and settings"), paste0(td, "/README.txt"))
        zip::zip(file, list.files(td, recursive = TRUE), root = td)
      })
    })

  session$onSessionEnded(function() { j <- isolate(rv$job); if (!is.null(j$job) && j$job$is_alive()) j$job$kill() })
}

shinyApp(ui, server)
