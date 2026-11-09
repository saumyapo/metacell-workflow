# =============================================================================
#  MCexplorer  |  global.R
#  Shiny loads this first: settings, packages, and the helpers used by ui.R and server.R.
#
#  START
#    RStudio:  open global.R, ui.R or server.R and click "Run App".
#    R:        shiny::runApp("path/to/MCexplorer")
#    Minerva:  MCX_DATA=/path/to/folder Rscript -e 'shiny::runApp("MCexplorer", port = 8888, launch.browser = FALSE)'
#
#  FILES
#    Needed (from the metacell gene-module step):
#      ...genelist...rd   gene modules               (object: genelist)
#      ...metadata...rd   module scores per metacell (object: metadf_z)
#      ...counts...rd     summed counts per metacell (object: dgcmtxcounts)
#      ...HMlist...rd     the module heatmap         (object: ht)
#    Optional:
#      cell information   annotated Seurat object (.rds) or a table of cells (.csv / .tsv) with a metacell column
#      metacell UMAP      .csv with metacell names and two coordinate columns
#
#  APP FILES
#    global.R    settings, packages, helpers        ui.R      layout
#    server.R    everything that reacts             help.R    text of the "How to use" windows
#    genesets.R  gene-set library (MSigDB via msigdbr, or .gmt files)
#    assets/     stylesheet, browser scripts and SortableJS (drag-and-drop lists)
#
#  OPTIONAL PACKAGES
#    callr (background work), uwot (UMAP), msigdbr (gene-set library)
# =============================================================================

DATA_PATH <- ""   # optional: folder that holds the files above
if (!nzchar(DATA_PATH)) DATA_PATH <- Sys.getenv("MCX_DATA", "")

local({
  need <- c(shiny = "1.7.5", bslib = "0.6.0", ggplot2 = "3.4.0", patchwork = "1.1.0", ComplexHeatmap = "2.14.0",
            circlize = "0.4.15", Matrix = "1.5.0", DT = "0.27", colourpicker = "1.2.0", openxlsx = "4.2.0")
  have <- vapply(names(need), function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, "")
  bad <- names(need)[is.na(have) | mapply(function(h, n) !is.na(h) && utils::compareVersion(h, n) < 0, have, need)]
  if (length(bad)) stop("MCexplorer needs these packages installed or updated: ",
                        paste0(bad, " (have ", ifelse(is.na(have[bad]), "none", have[bad]), ", need >= ", need[bad], ")", collapse = ", "),
                        "\n  install.packages(c(", paste0("\"", setdiff(bad, "ComplexHeatmap"), "\"", collapse = ", "), "))",
                        if ("ComplexHeatmap" %in% bad) "\n  BiocManager::install(\"ComplexHeatmap\")" else "", call. = FALSE)
})

for (p in c("callr", "uwot", "msigdbr")) if (!requireNamespace(p, quietly = TRUE))
  message("Optional package not installed: ", p, switch(p,
          callr = " (cell information and new blocks are computed in the background with it)", uwot = " (needed to compute a UMAP)",
          msigdbr = " (gene-set library: Hallmark, GO, Reactome, KEGG)"))
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
addResourcePath("assets", normalizePath("assets"))   # stylesheet and scripts of the app

APP_VERSION  <- "10.0"
EPS          <- 1e-5
CACHE_VER    <- 6
AUTOSAVE_DIR <- tools::R_user_dir("MCexplorer", "data")
ACCENT       <- "#0F766E"

PAL <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#86BCB6",
         "#F1CE63","#D37295","#A0CBE8","#FFBE7D","#8CD17D","#B6992D","#499894","#FABFD2","#D4A6C8","#2F4B7C",
         "#A05195","#665191","#D45087","#79706E")
PALETTES <- list(
  "Tableau"     = PAL,
  "Okabe-Ito"   = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999","#000000"),
  "Set2 (soft)" = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3"),
  "Dark2"       = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"),
  "Pastel"      = c("#AEC7E8","#FFBB78","#98DF8A","#FF9896","#C5B0D5","#C49C94","#F7B6D2","#DBDB8D","#9EDAE5","#C7C7C7"),
  "Viridis"     = "viridis",
  "Greys"       = "grays")
GREY <- "#C9CDD4"
NO_LABEL <- c("undefined", "Undefined", "Unassigned", "Outlier", "NA")

`%|na|%` <- function(a, b) if (length(a) && !is.na(a[1])) a else b
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ============================================================== small helpers

mcx_log  <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
fmt_secs <- function(s) { s <- as.numeric(s); if (s < 60) sprintf("%.0f s", s) else sprintf("%d min %02d s", as.integer(s %/% 60), as.integer(round(s %% 60))) }
fmt_size <- function(b) format(structure(b, class = "object_size"), units = "auto")
fmt_n    <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
pretty_name <- function(x) gsub("_", " ", x)
safe_name <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))

parse_list <- function(txt) {
  if (is.null(txt) || !nzchar(trimws(paste(txt, collapse = " ")))) return(character(0))
  s <- unlist(strsplit(gsub("^\\s*c\\(|\\)\\s*$", "", trimws(paste(txt, collapse = "\n"))), "[,;\n\t ]+"))
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
hex6 <- function(x) { m <- grDevices::col2rgb(x); grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255) }
valid_col <- function(x) !is.na(x) & nzchar(x) & vapply(x, function(z) !inherits(try(grDevices::col2rgb(z), silent = TRUE), "try-error"), TRUE)
col_or <- function(x, def) { x <- trimws(x %||% ""); if (length(x) == 1 && valid_col(x)) x else def }
palette_cols <- function(name, n) {
  p <- PALETTES[[name]] %||% PAL
  if (identical(p, "viridis")) return(grDevices::hcl.colors(n, "viridis"))
  if (identical(p, "grays")) return(grDevices::gray.colors(n, 0.15, 0.8))
  rep_len(p, n)
}
# distinct colours beyond the palette: golden-angle hue steps, alternating lightness
spread_col <- function(n) grDevices::hcl(h = (n * 137.508) %% 360, c = 58, l = c(62, 72, 52)[floor(n) %% 3 + 1])
file_key <- function(paths) paste(basename(paths), file.size(paths), as.numeric(file.mtime(paths)), collapse = "|")
chunks <- function(n, size = 1000) split(seq_len(n), ceiling(seq_len(n) / size))

# ============================================================== reading files

load_rd <- function(path, expected) {
  e <- new.env(); nm <- load(path, envir = e)
  if (expected %in% nm) return(get(expected, envir = e))
  if (length(nm) == 1) return(get(nm, envir = e))
  stop(basename(path), " does not contain '", expected, "'")
}
rd_role <- function(obj) {
  if (inherits(obj, c("HeatmapList", "Heatmap"))) return("heatmap")
  if (inherits(obj, "Matrix") || (is.matrix(obj) && is.numeric(obj))) return("counts")
  if (is.data.frame(obj) && any(grepl("_scaled$", colnames(obj)))) return("scores")
  if (is.data.frame(obj) && ncol(obj) >= 2 && is.character(obj[[2]])) return("genelist")
  NA_character_
}
ROLE_KEYS  <- c(genelist = "genelist", scores = "metadata", counts = "counts", heatmap = "hmlist")
ROLE_NAMES <- c(genelist = "Gene modules", scores = "Module scores", counts = "Counts per metacell", heatmap = "Module heatmap")

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
  if (!length(cand)) return(NULL)
  a <- cand[grepl("annot", basename(cand), ignore.case = TRUE)]
  if (!length(a)) return(cand[1])
  tok <- strsplit(stem, "_")[[1]]
  a[which.max(vapply(a, function(x) sum(vapply(tok, grepl, TRUE, x = tolower(basename(x)), fixed = TRUE)), 0))]
}
read_umap <- function(path, mcs) {
  u <- if (grepl("\\.csv$", path, ignore.case = TRUE)) utils::read.csv(path, check.names = FALSE) else utils::read.delim(path, check.names = FALSE)
  idc <- which(vapply(u, function(v) sum(as.character(v) %in% mcs), 0) > 0.5 * length(mcs))
  ids <- if (length(idc)) as.character(u[[idc[1]]]) else rownames(u)
  num <- setdiff(which(vapply(u, is.numeric, TRUE)), idc)
  if (length(num) < 2 || sum(ids %in% mcs) < 0.5 * length(mcs)) stop("Could not find metacell names and two coordinate columns in this file.")
  xy <- data.frame(x = u[[num[1]]], y = u[[num[2]]], row.names = ids)
  xy[intersect(rownames(xy), mcs), , drop = FALSE]
}

# ============================================================== dataset preparation

reporter <- function(n, title) {
  i <- 0; t0 <- Sys.time(); tl <- Sys.time()
  mcx_log("Opening ", title)
  function(label) {
    if (i > 0) mcx_log(sprintf("      done (%s)", fmt_secs(difftime(Sys.time(), tl, units = "secs"))))
    i <<- i + 1; tl <<- Sys.time()
    mcx_log(sprintf("[%d/%d] %s", i, n, label))
    if (!is.null(getDefaultReactiveDomain())) incProgress(1 / (n + 1), detail = label)
    invisible(t0)
  }
}

# per-group sums for top genes, one pass over the matrix
# S: sum of log2(1 + fraction * 1e5) (= log2(fraction + 1e-5) shifted); frac: mean fraction
group_stats <- function(counts, libsize, grp, ng) {
  n <- ncol(counts); g <- nrow(counts)
  B <- sparseMatrix(i = seq_len(n), j = grp, x = 1, dims = c(n, ng))
  S <- Fr <- matrix(0, g, ng); s1 <- s2 <- numeric(g)
  for (idx in chunks(n)) {
    m <- as(counts[, idx, drop = FALSE] %*% Diagonal(x = 1 / libsize[idx]), "CsparseMatrix")
    Bi <- B[idx, , drop = FALSE]
    Fr <- Fr + as.matrix(m %*% Bi)
    m@x <- log2(m@x + EPS) - log2(EPS)
    S <- S + as.matrix(m %*% Bi)
    s1 <- s1 + Matrix::rowSums(m); m@x <- m@x^2; s2 <- s2 + Matrix::rowSums(m)
  }
  sizes <- tabulate(grp, ng)
  dimnames(S) <- dimnames(Fr) <- list(rownames(counts), seq_len(ng))
  list(S = S, frac = sweep(Fr, 2, pmax(sizes, 1), "/"), sizes = sizes,
       gene_frac = rowSums(Fr) / n, gene_mean = setNames(s1 / n, rownames(counts)),
       gene_var = setNames((s2 - s1^2 / n) / max(n - 1, 1), rownames(counts)))
}
rank_genes <- function(gs, A, B = NULL, min_frac = EPS) {
  sz <- gs$sizes; n <- sum(sz)
  SA <- rowSums(gs$S[, A, drop = FALSE]); nA <- sum(sz[A])
  if (is.null(B) || !length(B)) { SB <- rowSums(gs$S) - SA; nB <- n - nA } else { SB <- rowSums(gs$S[, B, drop = FALSE]); nB <- sum(sz[B]) }
  fa <- as.vector(gs$frac[, A, drop = FALSE] %*% sz[A]) / max(nA, 1)
  out <- data.frame(gene = rownames(gs$S), log2FC = SA / max(nA, 1) - SB / max(nB, 1), mean_frac = fa, stringsAsFactors = FALSE)
  out <- out[out$mean_frac >= min_frac, ]
  out[order(-out$log2FC), ]
}
top_genes_all <- function(gs, n = 10, min_frac = 5e-5) {
  sz <- gs$sizes; tot <- rowSums(gs$S); N <- sum(sz)
  lfc <- sweep(gs$S, 2, pmax(sz, 1), "/") - sweep(tot - gs$S, 2, pmax(N - sz, 1), "/")
  lapply(seq_along(sz), function(k) { ok <- gs$frac[, k] >= min_frac; g <- rownames(gs$S)[ok]; head(g[order(-lfc[ok, k])], n) })
}

prepare_dataset <- function(src, name, cache_file = NULL) {
  step <- reporter(6, name)
  t0 <- step("Reading gene modules")
  genelist <- if (is.character(src$genelist)) load_rd(src$genelist, "genelist") else src$genelist
  step("Reading module scores")
  metadf_z <- if (is.character(src$scores)) load_rd(src$scores, "metadf_z") else src$scores
  step("Reading the saved heatmap")
  ht <- if (is.character(src$heatmap)) load_rd(src$heatmap, "ht") else src$heatmap
  step(paste0("Reading counts", if (is.character(src$counts)) paste0(" (", fmt_size(file.size(src$counts)), ", the largest file)") else ""))
  counts <- if (is.character(src$counts)) load_rd(src$counts, "dgcmtxcounts") else src$counts
  counts <- as(counts, "CsparseMatrix")

  step("Checking that the files match")
  genelist <- as.data.frame(genelist); colnames(genelist)[1:2] <- c("Mods", "Gens")
  zcols <- grep("_scaled$", colnames(metadf_z), value = TRUE)
  if (!length(zcols)) stop("The module scores file has no '<module>_scaled' columns.")
  Z <- as.matrix(metadf_z[, zcols, drop = FALSE])
  if (!identical(rownames(Z), colnames(counts))) {
    if (length(intersect(rownames(Z), colnames(counts))) < nrow(Z)) stop("Metacell names differ between the module scores and counts files.")
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
  glist <- lapply(strsplit(gtxt, ","), function(g) trimws(g[nzchar(trimws(g))]))
  mods <- data.frame(id = zcols, num = num, stringsAsFactors = FALSE)
  mods$name  <- paste("Module", num)
  mods$short <- paste0("M", num)
  mods$genes <- vapply(glist, paste, "", collapse = ", ")
  mods$n_genes <- lengths(glist)
  mods$lab   <- paste0(mods$short, "  ", vapply(glist, function(g) paste(head(g, 3), collapse = " "), ""))
  mods$title <- paste0(mods$name, ": ", vapply(glist, function(g) paste(head(g, 5), collapse = ", "), ""), ifelse(mods$n_genes > 5, ", ...", ""))
  libsize <- Matrix::colSums(counts)

  ckey <- if (is.character(src$counts)) file_key(unlist(src[c("genelist", "scores", "counts", "heatmap")])) else NULL
  cache <- cache_read(cache_file, ckey)
  if (!is.null(cache$bs)) { step("Summarising blocks (from cache)"); bs <- cache$bs }
  else { step("Summarising blocks (once; cached for next time)"); bs <- group_stats(counts, libsize, block, nb); cache_write(cache_file, ckey, "bs", bs) }
  mcx_log(sprintf("      ready in %s: %s metacells, %d modules, %d blocks", fmt_secs(difftime(Sys.time(), t0, units = "secs")), fmt_n(length(mcs)), ncol(Z), nb))

  list(name = name, src = lapply(src, function(x) if (is.character(x)) x else "(uploaded)"),
       Z = Z, counts = counts, ht = ht, row_ord = row_ord, col_ord = col_ord,
       mcs = mcs, nb = nb, block = block, mc_order_ht = mcs[unlist(col_ord)], mods = mods, mod_genes = setNames(glist, zcols),
       gene2mod = setNames(rep(mods$short, mods$n_genes), unlist(glist)),
       libsize = libsize, nfeat = diff(counts@p), bs = bs,
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

# New blocks, made exactly the way the pipeline makes its blocks: ComplexHeatmap's column_km, i.e. k-means of
# the metacells on their module z-scores (iter.max 50), blocks ordered by clustering of their means, and
# metacells ordered inside each block by complete-linkage clustering (Gu 2016, Bioinformatics 32:2847).
# The pipeline calls set.seed(42) right before Heatmap(), so the same seed and number of blocks gives the
# pipeline's blocks back exactly (checked against the pipeline call on 1,200 and 15,000 metacells).
# Every metacell stays whole: blocks only group metacells.
km_blocks <- function(Z, k, seed = 42) {
  k <- max(2, min(as.integer(k), nrow(Z) - 1))
  set.seed(seed)
  ht <- suppressMessages(ComplexHeatmap::prepare(ComplexHeatmap::Heatmap(t(Z), cluster_rows = TRUE, row_dend_reorder = FALSE,
          cluster_columns = TRUE, column_dend_reorder = FALSE, column_km = k, show_column_names = FALSE, show_row_names = FALSE)))
  co <- ComplexHeatmap::column_order(ht); if (!is.list(co)) co <- list(co)
  id <- integer(nrow(Z)); for (i in seq_along(co)) id[co[[i]]] <- i
  list(id = id, order = unlist(co, use.names = FALSE))
}

# Runs in a background R process: reads cell information and summarises it per metacell.
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
  labels <- list(); comp <- list(); comp_major <- list(); comp_purity <- list(); nums <- list(); lv <- list(); total_umis <- NULL; fac <- list()
  for (c in setdiff(names(md), c(mc_col, "..mc"))) {
    v <- md[[c]]
    if (is.numeric(v) && !is.factor(v)) {
      if (!cell_level || length(unique(v[seq_len(min(5000, length(v)))])) <= 1) next
      s <- as.vector(rowsum(as.numeric(v), mcf, reorder = TRUE, na.rm = TRUE)[as.character(mcs), 1]); s[is.na(s)] <- 0
      if (is.null(total_umis) && grepl("^(ncount_rna|total_umis?|umis?)$", c, ignore.case = TRUE)) total_umis <- s
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
    else { comp[[c]] <- tab; comp_major[[c]] <- major; comp_purity[[c]] <- round(mx / pmax(rs, 1), 3); fac[[c]] <- f }
  }
  pairs <- list()
  for (a in names(fac)) for (z in setdiff(names(fac), a)) { t <- table(fac[[a]], fac[[z]]); pairs[[a]][[z]] <- matrix(t, nrow = nrow(t), dimnames = dimnames(t)) }
  list(file = path, mc_col = mc_col, cell_level = cell_level, n_total = nrow(md),
       n_cells = if (cell_level) n_cells else NULL, total_umis = total_umis,
       labels = labels, comp = comp, comp_major = comp_major, comp_purity = comp_purity, nums = nums, levels = lv, pairs = pairs)
}

# ============================================================== expression and scores

# log2(1 + fraction * 1e5) for some genes (genes x metacells, sparse); equals log2(fraction + 1e-5) shifted
log_expr <- function(b, genes, mcs = NULL) {
  X <- if (is.null(mcs)) b$counts[genes, , drop = FALSE] else b$counts[genes, mcs, drop = FALSE]
  ls <- if (is.null(mcs)) b$libsize else b$libsize[mcs]
  X <- as(X %*% Diagonal(x = 1e5 / ls), "CsparseMatrix"); X@x <- log2(1 + X@x)
  X
}
gene_expr <- function(b, g, mcs = NULL) as.vector(log_expr(b, g, mcs))

# signature score per metacell
#   zmean: each gene z-scored across metacells, then averaged (every gene counts equally)
#   ctrl:  mean of the genes minus mean of expression-matched control genes (Tirosh et al. 2016; Seurat AddModuleScore)
score_signature <- function(b, genes, method = "zmean", nbin = 24, nctrl = 100) {
  X <- log_expr(b, genes); n <- ncol(X)
  if (method == "ctrl") {
    gm <- b$bs$gene_mean; gm <- gm[is.finite(gm) & b$bs$gene_frac > 0]
    set.seed(1); bins <- cut(rank(gm + stats::rnorm(length(gm)) / 1e30, ties.method = "first"), nbin, labels = FALSE); names(bins) <- names(gm)
    pool <- setdiff(names(gm), genes)
    ctrl <- unique(unlist(lapply(genes, function(g) { bg <- bins[g]; if (is.na(bg)) return(NULL); cand <- pool[bins[pool] == bg]; sample(cand, min(nctrl, length(cand))) })))
    s <- Matrix::colMeans(X) - Matrix::colMeans(log_expr(b, ctrl))
  } else {
    mu <- Matrix::rowMeans(X); v <- Matrix::rowMeans(X^2) - mu^2; sd <- sqrt(pmax(v, 0) * n / max(n - 1, 1))
    ok <- sd > 0; if (!any(ok)) return(setNames(rep(0, n), b$mcs))
    s <- (Matrix::colSums(Diagonal(x = 1 / sd[ok]) %*% X[ok, , drop = FALSE]) - sum(mu[ok] / sd[ok])) / sum(ok)
  }
  setNames(as.vector(s), b$mcs)
}

# ============================================================== statistics

sig_stars <- function(p) ifelse(is.na(p), "", as.character(cut(p, c(-Inf, 1e-4, 1e-3, 0.01, 0.05, Inf), c("****", "***", "**", "*", "ns"), right = FALSE)))
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), as.character(signif(p, 2))))

p_two <- function(a, b, test = "wilcox") {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  if (test == "t") tryCatch(stats::t.test(a, b)$p.value, error = function(e) NA_real_)
  else tryCatch(suppressWarnings(stats::wilcox.test(a, b)$p.value), error = function(e) NA_real_)
}
p_overall <- function(y, g, test = "wilcox") {
  ok <- is.finite(y) & !is.na(g); y <- y[ok]; g <- as.character(g[ok])
  cnt <- table(g); keep <- g %in% names(cnt)[cnt >= 2]; y <- y[keep]; g <- factor(g[keep])
  if (nlevels(g) < 2) return(NA_real_)
  if (nlevels(g) == 2) { l <- levels(g); return(p_two(y[g == l[1]], y[g == l[2]], test)) }
  if (test == "t") tryCatch(summary(stats::aov(y ~ g))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
  else tryCatch(stats::kruskal.test(y, g)$p.value, error = function(e) NA_real_)
}
# heatmap and dot plot statistics: every column against the reference column, in every row
tile_tests <- function(df, ref, test = "wilcox", adj = "holm", lab = "p", hidens = FALSE) {
  out <- list()
  for (f in levels(droplevels(factor(df$feature)))) {
    d <- df[df$feature == f, ]; lv <- levels(droplevels(factor(d$x)))
    for (g in setdiff(lv, ref)) out[[length(out) + 1]] <- data.frame(Row = f, Group = g, Reference = ref, p = p_two(d$y[d$x == g], d$y[d$x == ref], test))
  }
  if (!length(out)) return(NULL)
  st <- do.call(rbind, out); st$`Adjusted p` <- stats::p.adjust(st$p, if (adj == "none") "none" else adj); st$Significance <- sig_stars(st$`Adjusted p`)
  st$label <- if (identical(lab, "stars")) st$Significance else fmt_p(st$`Adjusted p`)
  if (isTRUE(hidens)) st$label[is.na(st$`Adjusted p`) | st$`Adjusted p` >= 0.05] <- ""
  st
}
make_pairs <- function(ok, cmp, ref = NULL, pairs = NULL) {
  out <- switch(cmp %||% "all",
    all  = if (length(ok) > 1) utils::combn(ok, 2, simplify = FALSE) else list(),
    ref  = if (isTRUE(ref %in% ok)) lapply(setdiff(ok, ref), function(g) c(ref, g)) else list(),
    pick = lapply(strsplit(pairs %||% character(0), " vs ", fixed = TRUE), identity))
  Filter(function(z) length(z) == 2 && all(z %in% ok), out)
}

# tests within each panel (facet) and, when split, within each x group; returns brackets ready to draw
run_stats <- function(df, opt) {
  panels <- if (!is.null(df$facet)) levels(droplevels(factor(df$facet))) else NA_character_
  split <- isTRUE(opt$split); xl <- levels(df$x); sl <- levels(df$fill)
  over <- list(); prs <- list()
  for (pn in panels) {
    d0 <- if (is.na(pn)) df else df[df$facet == pn, ]
    parts <- if (split) split(d0, d0$x, drop = TRUE) else list(all = d0)
    for (bn in names(parts)) {
      d <- parts[[bn]]; d <- d[is.finite(d$y), ]
      g <- factor(as.character(if (split) d$fill else d$x), levels = if (split) sl else xl)
      cnt <- table(g); ok <- names(cnt)[cnt >= 2]
      if (length(ok) < 2) next
      po <- p_overall(d$y, g, opt$test)
      over[[length(over) + 1]] <- data.frame(Panel = if (is.na(pn)) "" else pn, Within = if (split) bn else "",
        Test = if (length(ok) == 2) (if (opt$test == "t") "t-test" else "Wilcoxon") else (if (opt$test == "t") "ANOVA" else "Kruskal-Wallis"),
        Groups = length(ok), p = po, stringsAsFactors = FALSE)
      comps <- make_pairs(ok, opt$cmp, opt$ref, opt$pairs); if (!length(comps)) next
      pv <- vapply(comps, function(z) p_two(d$y[g == z[1]], d$y[g == z[2]], opt$test), 0)
      pa <- stats::p.adjust(pv, opt$adj %||% "holm")
      rng <- diff(range(d0$y, na.rm = TRUE)); if (!is.finite(rng) || rng == 0) rng <- 1
      top <- max(d$y, na.rm = TRUE)
      if (split) { ns <- length(sl); ctr <- match(bn, xl); off <- function(l) (match(l, sl) - (ns + 1) / 2) * 0.8 / ns
                   x1 <- ctr + off(vapply(comps, `[`, "", 1)); x2 <- ctr + off(vapply(comps, `[`, "", 2)) }
      else { x1 <- match(vapply(comps, `[`, "", 1), xl); x2 <- match(vapply(comps, `[`, "", 2), xl) }
      prs[[length(prs) + 1]] <- data.frame(Panel = if (is.na(pn)) "" else pn, Within = if (split) bn else "",
        group1 = vapply(comps, `[`, "", 1), group2 = vapply(comps, `[`, "", 2),
        n1 = vapply(comps, function(z) sum(g == z[1]), 1L), n2 = vapply(comps, function(z) sum(g == z[2]), 1L),
        p = pv, p_adj = pa, xmin = pmin(x1, x2), xmax = pmax(x1, x2),
        y = top + rng * (0.07 + 0.1 * (seq_along(comps) - 1)), tip = rng * 0.015, stringsAsFactors = FALSE)
    }
  }
  pw <- if (length(prs)) do.call(rbind, prs) else NULL
  if (!is.null(pw)) { pw$stars <- sig_stars(pw$p_adj); pw$label <- fmt_p(pw$p_adj) }
  list(pw = pw, overall = if (length(over)) do.call(rbind, over) else NULL, adj = opt$adj, test = opt$test)
}
add_brackets <- function(p, st, opt, facet_col = NULL) {
  pw <- st$pw; if (is.null(pw)) return(p)
  if (isTRUE(opt$hidens)) pw <- pw[!is.na(pw$p_adj) & pw$p_adj < 0.05, ]
  if (!nrow(pw)) return(p)
  if (!is.null(facet_col)) pw$facet <- factor(pw$Panel, levels = facet_col)
  pw$lab <- if (identical(opt$label, "stars")) pw$stars else pw$label
  fc <- pw[, intersect("facet", names(pw)), drop = FALSE]
  seg <- rbind(data.frame(fc, x = pw$xmin, xend = pw$xmin, y = pw$y - pw$tip, yend = pw$y),
               data.frame(fc, x = pw$xmin, xend = pw$xmax, y = pw$y, yend = pw$y),
               data.frame(fc, x = pw$xmax, xend = pw$xmax, y = pw$y - pw$tip, yend = pw$y))
  pw$xm <- (pw$xmin + pw$xmax) / 2
  p + geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = 0.35, colour = "#3A3F4B") +
    geom_text(data = pw, aes(x = xm, y = y, label = lab), inherit.aes = FALSE, vjust = -0.3, size = opt$size %||% 3.5, colour = "#1D2129") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))
}
stats_table <- function(st) {
  out <- list()
  if (!is.null(st$overall)) out[[1]] <- data.frame(Panel = st$overall$Panel, Within = st$overall$Within,
    Comparison = paste0("all ", st$overall$Groups, " groups (", st$overall$Test, ")"), n = "", p = st$overall$p, `Adjusted p` = NA, Significance = sig_stars(st$overall$p), check.names = FALSE)
  if (!is.null(st$pw)) out[[2]] <- data.frame(Panel = st$pw$Panel, Within = st$pw$Within, Comparison = paste(st$pw$group1, "vs", st$pw$group2),
    n = paste(st$pw$n1, st$pw$n2, sep = " / "), p = st$pw$p, `Adjusted p` = st$pw$p_adj, Significance = st$pw$stars, check.names = FALSE)
  if (!length(out)) return(NULL)
  d <- do.call(rbind, out)
  d[, c(if (any(nzchar(d$Panel))) "Panel", if (any(nzchar(d$Within))) "Within", "Comparison", "n", "p", "Adjusted p", "Significance")]
}

# ============================================================== plots

theme_mcx <- function(base = 13) {
  theme_minimal(base_size = base) + theme(
    plot.title = element_text(face = "bold", size = rel(1.05), colour = "#1D2129", margin = margin(b = 4)),
    plot.subtitle = element_text(colour = "#5B6170", size = rel(0.8), margin = margin(b = 8)),
    axis.title = element_text(colour = "#3A3F4B", size = rel(0.88)),
    axis.text = element_text(colour = "#3A3F4B", size = rel(0.85)),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#ECEEF2", linewidth = 0.4),
    axis.line.x = element_line(colour = "#9AA0AC", linewidth = 0.4),
    strip.text = element_text(face = "bold", size = rel(0.9), colour = "#1D2129"),
    legend.title = element_text(size = rel(0.82), face = "bold", colour = "#3A3F4B"),
    legend.text = element_text(size = rel(0.8)), legend.key.size = unit(0.9, "lines"),
    plot.margin = margin(10, 14, 8, 8))
}
x_text <- function(labels, angle = "auto") {
  ang <- if (identical(angle, "auto")) (if (max(nchar(labels), 0, na.rm = TRUE) > 7 || length(labels) > 10) 45 else 0) else as.numeric(angle)
  element_text(angle = ang, hjust = if (ang > 0) 1 else 0.5, vjust = if (ang == 90) 0.5 else 1)
}
lab_fun <- function(disp) function(v) unname(ifelse(is.na(disp[v]), v, disp[v]))

# box / violin / bar / strip, with optional points; split = dodge by a second variable
dist_plot <- function(df, type, pal, opt) {
  sp <- isTRUE(opt$split); dodge <- position_dodge(width = 0.8)
  pts <- if (sp) position_jitterdodge(jitter.width = 0.18, dodge.width = 0.8, seed = 1) else position_jitter(width = 0.18, height = 0, seed = 1)
  fa <- opt$fill_alpha %||% 0.35
  p <- ggplot(df, aes(x, y, fill = fill, colour = fill))
  if (type == "box") p <- p + geom_boxplot(width = 0.7, alpha = fa, outlier.shape = NA, linewidth = 0.5, position = dodge)
  if (type == "violin") p <- p + geom_violin(scale = "width", width = 0.8, alpha = fa, linewidth = 0.4, position = dodge) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.4, position = dodge, show.legend = FALSE)
  if (type == "bar") p <- p + stat_summary(fun = mean, geom = "col", width = 0.7, alpha = fa + 0.2, position = dodge) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25, linewidth = 0.5, position = dodge)
  if (type == "strip") p <- p + stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5, linewidth = 0.35, position = dodge)
  if (isTRUE(opt$points) || type == "strip") p <- p + geom_point(position = pts, size = opt$pt %||% 1.2, alpha = opt$alpha %||% 0.6, shape = 16)
  if (isTRUE(opt$mean)) p <- p + geom_hline(yintercept = mean(df$y, na.rm = TRUE), linetype = 2, colour = "#9AA0AC", linewidth = 0.4)
  xl <- levels(df$x); disp <- opt$xlabels %||% setNames(xl, xl)
  p + scale_fill_manual(values = pal, drop = FALSE, labels = lab_fun(opt$fill_labels %||% disp)) +
    scale_colour_manual(values = pal, drop = FALSE, labels = lab_fun(opt$fill_labels %||% disp)) +
    scale_x_discrete(drop = FALSE, labels = lab_fun(disp)) + theme_mcx(opt$txt %||% 13) +
    theme(axis.text.x = x_text(disp[xl], opt$angle %||% "auto"), legend.position = if (isTRUE(opt$legend)) "right" else "none")
}

# heatmap or dot plot of means: features (rows) x groups (columns)
#   opt$size_by: "pos" = % of points above zero (from df$pos), "mean" = the mean itself, "above" = % above the feature's overall median
summary_plot <- function(df, type, opt) {
  f <- droplevels(factor(df$feature)); x <- df$x
  m <- tapply(df$y, list(f, x), mean, na.rm = TRUE); m <- m[, levels(x), drop = FALSE]
  sz <- switch(opt$size_by %||% "pos",
    pos = 100 * tapply(as.numeric(df$pos), list(f, x), mean, na.rm = TRUE),
    mean = m,
    above = { med <- tapply(df$y, f, stats::median, na.rm = TRUE); 100 * tapply(as.numeric(df$y > med[as.character(f)]), list(f, x), mean, na.rm = TRUE) })
  sz <- sz[, levels(x), drop = FALSE]
  scaled <- isTRUE(opt$scale) && ncol(m) > 1
  val <- if (scaled) { s <- t(scale(t(m))); s[!is.finite(s)] <- 0; s } else m
  rl <- rownames(m)
  d <- data.frame(feature = factor(rep(rl, ncol(m)), levels = rev(rl)), x = factor(rep(colnames(m), each = nrow(m)), levels = colnames(m)),
                  val = as.vector(val), raw = as.vector(m), size = as.vector(sz))
  disp <- opt$xlabels %||% setNames(colnames(m), colnames(m))
  p <- ggplot(d, aes(x, feature))
  lim <- if (scaled) max(abs(d$val), 0.5, na.rm = TRUE) else NULL
  fill_name <- if (scaled) "Scaled mean\n(z across\ncolumns)" else (opt$val_name %||% "Mean")
  if (type == "heatmap") {
    p <- p + geom_tile(aes(fill = val), colour = "white", linewidth = 0.5)
    if (isTRUE(opt$values)) p <- p + geom_text(aes(label = signif(raw, 2)), size = 3, colour = "#1D2129")
    p <- p + if (scaled) scale_fill_gradient2(low = opt$low %||% "#000080", mid = opt$mid %||% "#FFFFFF", high = opt$high %||% "#B22222", limits = c(-lim, lim), name = fill_name)
             else scale_fill_viridis_c(name = fill_name)
  } else {
    p <- p + geom_point(aes(size = size, colour = val)) + scale_size_area(max_size = opt$max_size %||% 9, name = opt$size_name %||% "%")
    p <- p + if (scaled) scale_colour_gradient2(low = opt$low %||% "#000080", mid = "#E9EAEE", high = opt$high %||% "#B22222", limits = c(-lim, lim), name = fill_name)
             else scale_colour_viridis_c(name = fill_name)
  }
  tl <- opt$tile_labels   # data.frame(Row, Group, label): result of each column against the reference
  if (!is.null(tl) && nrow(tl)) {
    tl <- tl[nzchar(tl$label) & tl$Group %in% colnames(m) & tl$Row %in% rl, , drop = FALSE]
    if (nrow(tl)) p <- p + geom_text(data = data.frame(x = factor(tl$Group, levels = colnames(m)), feature = factor(tl$Row, levels = rev(rl)), label = tl$label),
                                     aes(x, feature, label = label), inherit.aes = FALSE, size = (opt$statsize %||% 3.5) * 0.9,
                                     colour = "#111827", fontface = "bold", nudge_y = if (type == "dotplot") 0.36 else 0)
  }
  if (!is.null(opt$ref)) disp[opt$ref] <- paste0(disp[opt$ref], "\n(reference)")
  p + scale_x_discrete(labels = lab_fun(disp)) + labs(x = NULL, y = NULL, title = opt$title) +
    theme_mcx(opt$txt %||% 13) + theme(panel.grid.major.y = element_blank(), axis.line.x = element_blank(),
      axis.text.x = x_text(disp[colnames(m)], opt$angle %||% "auto"), legend.position = "right", legend.box = "horizontal",
      legend.box.just = "top") +
    guides(fill = guide_colourbar(barheight = unit(6, "lines")), colour = guide_colourbar(barheight = unit(6, "lines")))
}

mat_colors <- function(M, cols, lim) {
  pal <- grDevices::colorRampPalette(cols)(256)
  i <- round((pmin(pmax(M, lim[1]), lim[2]) - lim[1]) / (lim[2] - lim[1]) * 255) + 1; i[is.na(i)] <- 1
  matrix(pal[i], nrow(M))
}
legend_layer <- function(lim) geom_point(data = data.frame(x = 1, y = 1, v = lim), aes(x, y, fill = v), shape = 21, size = 0, stroke = 0, inherit.aes = FALSE)

# Module heatmap as one full-bleed ggplot: the drawing area is the whole image, so the browser can turn a mouse
# position into a column itself (assets/app.js): hover, click and drag never depend on Shiny's plot coordinates.
HM_CPX <- 5.4    # width of one label character in px (7.5 pt)
heat_layout <- function(rowlab, grp, nbars, W, H, labels = TRUE) {
  nc <- length(grp); W <- max(300, W); H <- max(200, H)
  left <- max(nchar(rowlab), 3) * HM_CPX + 12; right <- 92; ux <- nc / max(50, W - left - right)
  r <- rle(as.character(grp)); ends <- cumsum(r$lengths); mids <- ends - r$lengths / 2 + 0.5; colpx <- 1 / ux
  vert <- sum(r$lengths * colpx >= 11) > 14
  wide <- if (vert) r$lengths * colpx >= 11 else r$lengths * colpx >= nchar(r$values) * HM_CPX + 6
  show_lab <- isTRUE(labels) && length(r$values) > 1 && sum(wide) <= 90
  top <- if (!show_lab) 6 else if (vert) max(nchar(r$values[wide])) * HM_CPX + 10 else 18
  bh <- 11; bg <- 3; y0 <- 4; y1 <- H - top - nbars * (bh + bg) - 2
  list(W = W, H = H, nc = nc, x0 = 0.5 - left * ux, x1 = nc + 0.5 + right * ux, ux = ux, r = r, ends = ends, mids = mids, vert = vert, wide = wide,
       show_lab = show_lab, top = top, bh = bh, bg = bg, y0 = y0, y1 = y1, ytop = y1 + nbars * (bh + bg))
}
heat_plot <- function(M, grp, bars, bar_cols, opt) {
  nm <- nrow(M); nc <- ncol(M); nb <- length(bars)
  L <- heat_layout(c(rownames(M), names(bars)), grp, nb, opt$width_px %||% 900, opt$height %||% 500, opt$labels)
  ux <- L$ux; tpt <- 7.5; y0 <- L$y0; y1 <- L$y1; rowh <- (y1 - y0) / nm
  # more metacells than screen pixels: average neighbouring columns per pixel (bars: most common value) instead of
  # letting the image drop columns, so narrow blocks are drawn from all their metacells
  npx <- floor(nc / L$ux * (opt$pr %||% 2))   # device pixels across the heatmap
  Mr <- M; barsr <- bars
  if (nc > npx && npx > 10) {
    bin <- ceiling(seq_len(nc) * npx / nc)
    Mr <- t(rowsum(t(M), bin, reorder = TRUE) / tabulate(bin, npx))
    barsr <- lapply(bars, function(v) { v <- as.character(v); v[is.na(v)] <- "NA"
      vapply(split(v, bin), function(z) { t <- table(z); names(t)[which.max(t)] }, "") })
  }
  p <- ggplot() + annotation_raster(grDevices::as.raster(mat_colors(Mr, c(opt$low, opt$mid, opt$high), c(-opt$lim, opt$lim))),
                                    xmin = 0.5, xmax = nc + 0.5, ymin = y0, ymax = y1, interpolate = FALSE)
  if (rowh >= 3) p <- p + geom_text(data = data.frame(y = y1 - (seq_len(nm) - 0.5) * rowh, l = rownames(M)), aes(y = y, label = l), x = 0.5 - 5 * ux,
                                    hjust = 1, size = max(4, min(8, rowh * 0.64)) / .pt, colour = "#3A3F4B", inherit.aes = FALSE)
  lab_df <- data.frame(y = numeric(0), l = character(0)); yb <- y1 + L$bg
  for (j in seq_len(nb)) {
    v <- as.character(barsr[[j]]); v[is.na(v)] <- "NA"; cl <- unname(bar_cols[[j]][v]); cl[is.na(cl)] <- GREY
    p <- p + annotation_raster(grDevices::as.raster(matrix(cl, 1)), xmin = 0.5, xmax = nc + 0.5, ymin = yb, ymax = yb + L$bh, interpolate = FALSE)
    lab_df <- rbind(lab_df, data.frame(y = yb + L$bh / 2, l = names(bars)[j])); yb <- yb + L$bh + L$bg
  }
  ends <- L$ends
  if (isTRUE(opt$sep) && length(ends) > 1 && length(ends) <= 400 && nc > length(ends))
    p <- p + geom_segment(data = data.frame(x = ends[-length(ends)] + 0.5), aes(x = x, xend = x), y = y0, yend = yb - L$bg, colour = "#1D2129", linewidth = 0.25, inherit.aes = FALSE)
  if (nrow(lab_df)) p <- p + geom_text(data = lab_df, aes(y = y, label = l), x = 0.5 - 5 * ux, hjust = 1, size = tpt / .pt, colour = "#3A3F4B", inherit.aes = FALSE)
  if (L$show_lab) p <- p + geom_text(data = data.frame(x = L$mids[L$wide], l = L$r$values[L$wide]), aes(x, label = l), y = yb + 2, size = tpt / .pt,
                                     angle = if (L$vert) 90 else 0, hjust = if (L$vert) 0 else 0.5, vjust = if (L$vert) 0.5 else 0, colour = "#3A3F4B", inherit.aes = FALSE)
  lh <- min(130, (y1 - y0) * 0.8); ym <- (y0 + y1) / 2; lx <- nc + 0.5 + 16 * ux
  key <- matrix(mat_colors(matrix(seq(opt$lim, -opt$lim, length.out = 60)), c(opt$low, opt$mid, opt$high), c(-opt$lim, opt$lim)), ncol = 1)
  p <- p + annotation_raster(grDevices::as.raster(key), xmin = lx, xmax = lx + 12 * ux, ymin = ym - lh / 2, ymax = ym + lh / 2, interpolate = TRUE) +
    geom_text(data = data.frame(y = ym + c(-1, 0, 1) * lh / 2, l = c(-opt$lim, 0, opt$lim)), aes(y = y, label = l), x = lx + 16 * ux, hjust = 0, size = tpt / .pt, colour = "#3A3F4B", inherit.aes = FALSE) +
    annotate("text", x = lx, y = ym + lh / 2 + 8, label = "Module\nscore (z)", hjust = 0, vjust = 0, size = 8 / .pt, fontface = "bold", colour = "#1F2937")
  p + coord_cartesian(xlim = c(L$x0, L$x1), ylim = c(0, L$H), expand = FALSE, clip = "off") + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}
heat_panel <- function(M, cols, lim, legend, ysize, ends = NULL, sep = TRUE, top_labels = NULL) {
  nm <- nrow(M); nc <- ncol(M)
  p <- ggplot() + geom_blank(data = data.frame(x = c(0.5, nc + 0.5), y = c(0.5, nm + 0.5)), aes(x, y)) +
    annotation_raster(grDevices::as.raster(mat_colors(M, cols, lim)), xmin = 0.5, xmax = nc + 0.5, ymin = 0.5, ymax = nm + 0.5, interpolate = FALSE) +
    legend_layer(lim) + scale_fill_gradientn(colours = cols, limits = lim, name = legend)
  if (sep && length(ends) > 1 && length(ends) <= 400 && nc > length(ends)) p <- p + geom_vline(xintercept = ends[-length(ends)] + 0.5, colour = "#1D2129", linewidth = 0.25)
  p + scale_x_continuous(limits = c(0.5, nc + 0.5), expand = c(0, 0), position = "top", breaks = top_labels$at %||% NULL, labels = top_labels$lab %||% NULL) +
    scale_y_continuous(limits = c(0.5, nm + 0.5), expand = c(0, 0), breaks = nm:1, labels = rownames(M)) +
    theme_minimal(base_size = 11) +
    theme(axis.title = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_text(size = ysize, colour = "#3A3F4B"),
          axis.text.x.top = element_text(size = 8, colour = "#3A3F4B", angle = if (length(top_labels$at) > 16) 90 else 0, hjust = if (length(top_labels$at) > 16) 0 else 0.5, vjust = 0.5),
          legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8), plot.margin = margin(2, 10, 2, 6))
}
bar_panel <- function(bars, bar_cols, ysize) {
  nc <- length(bars[[1]]); nbar <- length(bars)
  img <- do.call(rbind, lapply(seq_along(bars), function(j) { v <- as.character(bars[[j]]); v[is.na(v)] <- "NA"; cl <- unname(bar_cols[[j]][v]); cl[is.na(cl)] <- GREY; cl }))
  ggplot() + geom_blank(data = data.frame(x = c(0.5, nc + 0.5), y = c(0.5, nbar + 0.5)), aes(x, y)) +
    annotation_raster(grDevices::as.raster(img), xmin = 0.5, xmax = nc + 0.5, ymin = 0.5, ymax = nbar + 0.5, interpolate = FALSE) +
    scale_x_continuous(limits = c(0.5, nc + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, nbar + 0.5), expand = c(0, 0), breaks = nbar:1, labels = names(bars)) +
    theme_minimal(base_size = 11) + theme(axis.title = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_blank(), axis.text.y = element_text(size = ysize + 1, colour = "#3A3F4B"), plot.margin = margin(6, 10, 0, 6))
}

save_plot <- function(file, draw, w, h, fmt) {
  w <- max(as.numeric(w %||% 8), 1); h <- max(as.numeric(h %||% 6), 1)
  if (fmt == "png") grDevices::png(file, width = w, height = h, units = "in", res = 300, type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
  else if (fmt == "svg") { if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite(file, width = w, height = h) else grDevices::svg(file, width = w, height = h) }
  else grDevices::pdf(file, width = w, height = h, useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  draw()
}

# ============================================================== Excel

write_xlsx_sheets <- function(file, sheets) {
  wb <- openxlsx::createWorkbook()
  hs <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#EEF2FF", border = "bottom", borderColour = "#C7CBD4")
  for (nm in names(sheets)) {
    d <- sheets[[nm]]
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, d, headerStyle = hs)
    openxlsx::freezePane(wb, nm, firstRow = TRUE)
    w <- pmin(60, pmax(10, nchar(names(d)) + 2, vapply(d, function(x) max(c(0, nchar(as.character(x))), na.rm = TRUE), 1) + 2))
    openxlsx::setColWidths(wb, nm, cols = seq_along(d), widths = w)
    if ("Colour" %in% names(d)) for (i in seq_len(nrow(d)))
      openxlsx::addStyle(wb, nm, openxlsx::createStyle(fgFill = d$Colour[i], fontColour = d$Colour[i]), rows = i + 1, cols = match("Colour", names(d)))
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}
read_annotation_file <- function(path, mcs) {
  x <- if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    sh <- openxlsx::getSheetNames(path); openxlsx::read.xlsx(path, sheet = if ("Metacells" %in% sh) "Metacells" else 1, check.names = FALSE)
  } else utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  idc <- which.max(vapply(x, function(v) sum(as.character(v) %in% mcs), 0))
  ac <- intersect(c("annotation", "Annotation", "new_annotation"), names(x))[1]
  if (!length(idc) || sum(as.character(x[[idc]]) %in% mcs) == 0 || is.na(ac)) stop("The file needs a column of metacell names and a column called 'annotation'.")
  data.frame(metacell = as.character(x[[idc]]), annotation = as.character(x[[ac]]), stringsAsFactors = FALSE)
}

source("help.R", local = TRUE)        # text of the "How to use" windows
source("genesets.R", local = TRUE)    # gene-set library

# ============================================================== UI building blocks

LOGO <- HTML('<svg viewBox="0 0 18 18" aria-hidden="true"><rect x="0" y="0" width="5" height="5" rx="1" fill="#0F766E"/><rect x="6.5" y="0" width="5" height="5" rx="1" fill="#8FCBC1"/><rect x="13" y="0" width="5" height="5" rx="1" fill="#D5DAE1"/><rect x="0" y="6.5" width="5" height="5" rx="1" fill="#D5DAE1"/><rect x="6.5" y="6.5" width="5" height="5" rx="1" fill="#0F766E"/><rect x="13" y="6.5" width="5" height="5" rx="1" fill="#E3A33B"/><rect x="0" y="13" width="5" height="5" rx="1" fill="#8FCBC1"/><rect x="6.5" y="13" width="5" height="5" rx="1" fill="#D5DAE1"/><rect x="13" y="13" width="5" height="5" rx="1" fill="#0F766E"/></svg>')

RM_BTN <- list(plugins = list("remove_button"))   # duplicate x buttons from newer Shiny are hidden in styles.css
TIP_OPTS <- list(container = "body", boundary = "window", trigger = "hover")   # hover only: a clicked button does not keep its tooltip   # tooltips float above everything, never clipped by a panel
tip <- function(text) tooltip(span(class = "mcx-tip", icon("circle-question")), text, placement = "right", options = TIP_OPTS)
lab <- function(text, help = NULL) if (is.null(help)) text else span(text, tip(help))
section <- function(title, ..., help = NULL, class = NULL) div(class = paste("mcx-section", class), div(class = "mcx-section-title", title, if (!is.null(help)) tip(help)), ...)
hint <- function(...) div(class = "mcx-hint", ...)
with_tip <- function(el, text) tooltip(el, text, placement = "bottom", options = TIP_OPTS)
# number box with - and + buttons (handled in the browser), easier to set than a slider
num_step <- function(id, label, value, min = NA, max = NA, step = 1, help = NULL) div(class = "mcx-numstep",
  if (!is.null(label)) tags$label(class = "control-label", `for` = id, lab(label, help)),
  div(class = "mcx-stepper",
    tags$button(type = "button", class = "btn btn-light mcx-step", `data-for` = id, `data-dir` = "-1", icon("minus")),
    numericInput(id, NULL, value, min = min, max = max, step = step),
    tags$button(type = "button", class = "btn btn-light mcx-step", `data-for` = id, `data-dir` = "1", icon("plus"))))
# a collapsible card: title row always visible, content folds away
fold <- function(title, ..., icon_name = NULL, open = TRUE, tools = NULL, value = "x", id = NULL)
  accordion(id = id, open = if (open) value else FALSE, class = "mcx-fold",
            accordion_panel(title = div(class = "mcx-foldhead", span(title), if (!is.null(tools)) div(class = "mcx-foldtools", onclick = "event.stopPropagation()", tools)),
                            value = value, icon = if (!is.null(icon_name)) icon(icon_name), ...))
goto_btn <- function(id, text, ic = "arrow-right", title = NULL) actionButton(id, text, icon = icon(ic), class = "btn-link mcx-goto", title = title)

# radio group with own markup (same Shiny radio binding), so it looks the same with any bslib version
radio_markup <- function(id, choices, selected, cls, item) {
  if (is.null(names(choices))) names(choices) <- choices
  selected <- selected %||% unname(choices[[1]])
  div(id = id, class = paste("shiny-input-radiogroup shiny-input-container", cls), role = "radiogroup",
      div(class = "shiny-options-group", lapply(seq_along(choices), function(i) item(unname(choices[[i]]), names(choices)[i], identical(unname(choices[[i]]), selected)))))
}
# segmented buttons; tips: a short hint per option, shown on hover
seg_input <- function(id, choices, selected = NULL, label = NULL, help = NULL, tips = NULL) tagList(
  if (!is.null(label)) tags$label(class = "control-label", lab(label, help)),
  radio_markup(id, choices, selected, "mcx-seg", function(v, l, on)
    tags$label(class = "mcx-seg-opt", `data-tip` = tips[[v]], `data-base-tip` = tips[[v]],
               tags$input(type = "radio", name = id, value = v, checked = if (on) NA), tags$span(l))))
radio_tips <- function(id, label, choices, tips = list(), selected = NULL, help = NULL) tagList(
  tags$label(class = "control-label", lab(label, help)),
  radio_markup(id, choices, selected, "mcx-radiolist", function(v, l, on)
    div(class = "mcx-radio-row", tags$label(tags$input(type = "radio", name = id, value = v, checked = if (on) NA), tags$span(l)),
        if (!is.null(tips[[v]])) tip(tips[[v]]))))
# large option cards: [list(value, icon, title, sub), ...]
mode_tabs <- function(id, options, selected = NULL, cls = NULL) {
  ch <- setNames(vapply(options, `[[`, "", "value"), vapply(options, `[[`, "", "title"))
  info <- setNames(options, ch)
  radio_markup(id, ch, selected, paste("mcx-modes", cls), function(v, l, on) {
    o <- info[[v]]
    tags$label(class = "mcx-mode", tags$input(type = "radio", name = id, value = v, checked = if (on) NA),
      div(class = "mcx-mode-box", div(class = "mcx-mode-icon", icon(o$icon)), div(div(class = "mcx-mode-title", o$title), div(class = "mcx-mode-sub", o$sub))))
  })
}
# draggable, tickable list (value: ticked items in display order); items are sent from the server
rank_input <- function(id, label = NULL, help = NULL) div(class = "mcx-rank-wrap",
  if (!is.null(label)) tags$label(class = "control-label", lab(label, help)),
  div(id = id, class = "mcx-rank",
      div(class = "mcx-rank-head", span(class = "mcx-rank-count"), span(class = "ms-auto"),
          tags$a(href = "#", class = "mcx-rank-all", "All"), tags$a(href = "#", class = "mcx-rank-none", "None")),
      div(class = "mcx-rank-list")))
# "All / None" links for a checkbox group
check_all <- function(id) div(class = "mcx-checkall",
  tags$a(href = "#", onclick = sprintf("$('#%s input[type=checkbox]').prop('checked', true).first().trigger('change'); return false;", id), "Select all"),
  tags$a(href = "#", onclick = sprintf("$('#%s input[type=checkbox]').prop('checked', false).first().trigger('change'); return false;", id), "None"))

# colour box (hex text + picker) with a colour-wheel button and an eyedropper
colour_field <- function(id, label, value) div(class = "mcx-colourfield",
  if (!is.null(label)) tags$label(class = "control-label", label),
  div(class = "mcx-colourrow",
    div(class = "mcx-colourbox", colourpicker::colourInput(id, NULL, value, showColour = "both")),
    tags$label(class = "btn btn-light mcx-iconbtn", title = "Colour wheel", icon("circle-half-stroke"),
               tags$input(type = "color", class = "mcx-native", `data-for` = id, value = hex6(value))),
    tags$button(type = "button", class = "btn btn-light mcx-iconbtn mcx-dropper", `data-for` = id, title = "Pick a colour from the screen", icon("eye-dropper"))))

# signature chooser: select + new (own gene list) + library (published gene sets) + edit
sig_picker <- function(id, multiple = FALSE) div(class = "mcx-sigpick",
  div(class = "mcx-sigpick-row",
    div(selectizeInput(id, NULL, choices = NULL, multiple = multiple,
        options = c(list(placeholder = "No signature yet"), if (multiple) RM_BTN))),
    with_tip(actionButton(paste0(id, "_new"), NULL, icon = icon("plus"), class = "btn-light"), "New signature from a typed or pasted gene list"),
    with_tip(actionButton(paste0(id, "_lib"), NULL, icon = icon("book-open"), class = "btn-light"), "Gene-set library: Hallmark, GO biological process, Reactome, KEGG, or a .gmt file"),
    with_tip(actionButton(paste0(id, "_edit"), NULL, icon = icon("pen"), class = "btn-light"), "Edit or delete the chosen signature")),
  uiOutput(paste0(id, "_info")))

# what to plot: module(s) / gene(s) / signature(s) / metric
FEATURE_TIPS <- list(module = "Module score: the average z-score of the module's genes in each metacell.",
                     gene = "Expression of single genes: log2 of the fraction of the metacell's counts.",
                     sig = "A gene list scored in each metacell. Made from a typed list or from the gene-set library.",
                     metric = "Properties of each metacell, such as its number of cells or UMIs.")
feature_picker <- function(p, metric = TRUE, multi = TRUE) tagList(
  seg_input(paste0(p, "_what"), c("Module" = "module", "Gene" = "gene", "Signature" = "sig", if (metric) c("Metric" = "metric")), tips = FEATURE_TIPS),
  conditionalPanel(sprintf("input.%s_what == 'module'", p),
    selectizeInput(paste0(p, "_mod"), NULL, choices = NULL, multiple = multi,
                   options = c(list(placeholder = "Search a module number or a gene"), if (multi) c(RM_BTN, list(maxItems = 20)))),
    uiOutput(paste0(p, "_modinfo"))),
  conditionalPanel(sprintf("input.%s_what == 'gene'", p),
    selectizeInput(paste0(p, "_gene"), NULL, choices = NULL, multiple = multi,
                   options = c(list(placeholder = if (multi) "Type one or more genes" else "Type a gene"), if (multi) c(RM_BTN, list(maxItems = 20))))),
  conditionalPanel(sprintf("input.%s_what == 'sig'", p), sig_picker(paste0(p, "_sig"), multiple = multi)),
  if (metric) conditionalPanel(sprintf("input.%s_what == 'metric'", p), selectInput(paste0(p, "_metric"), NULL, choices = NULL)),
  if (multi) hint("Several items are drawn as separate panels, or as rows of a heatmap or dot plot."))

SIG_METHODS <- c("Mean z-score of its genes" = "zmean", "Against matched control genes" = "ctrl")
SIG_METHOD_TIP <- paste("Mean z-score: each gene is scaled across metacells, then the genes are averaged (simple, every gene counts equally).",
                        "Matched controls: the average minus that of control genes with similar overall expression (Tirosh 2016; Seurat AddModuleScore), which removes differences in overall detection.")
PLOT_TYPES <- c("Box" = "box", "Violin" = "violin", "Bar" = "bar", "Strip" = "strip", "Heatmap" = "heatmap", "Dot plot" = "dotplot")
TYPE_TIPS <- list(box = "Median, quartiles and range", violin = "Shape of the distribution", bar = "Mean with standard error",
                  strip = "Every point with a median line", stacked = "Share of each group, stacked to 100%",
                  heatmap = "Mean per group; one row per item", dotplot = "Colour: mean; size: share of metacells above zero")

# plot type, points and statistics switches, shown in the plot card header
plot_toolbar <- function(p, types = PLOT_TYPES, stats_on = FALSE, points_on = TRUE, extra = NULL) div(class = "mcx-tools",
  seg_input(paste0(p, "_type"), types, tips = TYPE_TIPS),
  conditionalPanel(sprintf("['box','violin','bar'].indexOf(input.%s_type) > -1", p),
    with_tip(div(input_switch(paste0(p, "_points"), "Points", points_on)), "Show or hide the individual points.")),
  conditionalPanel(sprintf("['heatmap','dotplot'].indexOf(input.%s_type) > -1", p),
    div(id = paste0(p, "_scale_wrap"), with_tip(div(input_switch(paste0(p, "_scale"), "Scale rows", TRUE)), "Each row scaled across columns (z-score), so rows with different ranges can be compared."))),
  with_tip(div(input_switch(paste0(p, "_stats"), "Statistics", stats_on)), "Adds p-values to the plot and a table under it. The settings open in the left panel."),
  extra,
  with_tip(actionButton(paste0(p, "_colors"), NULL, icon = icon("palette"), class = "btn-light btn-sm mcx-colbtn"), "Colours of the groups"))

STAT_TIPS <- list(
  cmp = list(all = "Every group against every other. Gets crowded past 3 or 4 groups.",
             ref = "Each group against one reference group, e.g. every condition against a control.",
             pick = "Only the pairs added below."),
  test = list(wilcox = "Rank-based (Mann-Whitney U). No assumption about the shape of the values; the usual choice for few samples.",
              t = "Welch's t-test. Assumes roughly normal values; more power when that holds."),
  adj = list(holm = "Controls the chance of any false positive, less strict than Bonferroni. Good default.",
             BH = "Benjamini-Hochberg: controls the share of false positives among the significant results (FDR). Best when many comparisons are made.",
             bonferroni = "Multiplies each p-value by the number of comparisons. The strictest.",
             none = "Raw p-values. Only for a single planned comparison."))
# statistics settings: a normal collapsible panel, shown only while 'Statistics' is switched on
opt_hints <- function(p, field, tips) lapply(names(tips), function(v)
  conditionalPanel(sprintf("input.%s_%s == '%s'", p, field, v), div(class = "mcx-hint", tips[[v]])))
# hover text on a dropdown's label: every option with what it does
tip_options <- function(choices, tips) tip(HTML(paste0("<b>", names(choices), "</b>: ", vapply(unname(choices), function(v) tips[[v]], ""), collapse = "<br>")))
CMP_CH  <- c("All pairs" = "all", "Each against a reference" = "ref", "Chosen pairs" = "pick")
TEST_CH <- c("Wilcoxon" = "wilcox", "t-test" = "t")
ADJ_CH  <- c("Holm" = "holm", "Benjamini-Hochberg (FDR)" = "BH", "Bonferroni" = "bonferroni", "None" = "none")
stats_panel <- function(p, hide = FALSE) {
  heat <- sprintf("['heatmap','dotplot'].indexOf(input.%s_type) > -1", p)
  accordion_panel("Statistics settings", value = "stats", icon = icon("square-root-variable"),
    div(class = "mcx-hint", "Which groups are compared, with which test, and how p-values are corrected for testing many comparisons. Hover the ? next to each box for all its options."),
    conditionalPanel(paste0("!(", heat, ")"),
      tags$label(class = "control-label", "Comparisons", tip_options(CMP_CH, STAT_TIPS$cmp)),
      selectInput(paste0(p, "_cmp"), NULL, CMP_CH), opt_hints(p, "cmp", STAT_TIPS$cmp)),
    conditionalPanel(heat, div(class = "mcx-why", icon("circle-info"),
      span("Heatmap and dot plot: in every row, each column is tested against the reference column, and the result is written on the tile."))),
    conditionalPanel(sprintf("input.%s_cmp == 'ref' || %s", p, heat),
      selectInput(paste0(p, "_ref"), lab("Reference group", "The group every other group is compared with, e.g. the control condition."), choices = NULL)),
    conditionalPanel(sprintf("input.%s_cmp == 'pick' && !(%s)", p, heat),
      div(class = "mcx-pairpick",
        selectInput(paste0(p, "_pa"), "First group", choices = NULL), selectInput(paste0(p, "_pb"), "Compared with", choices = NULL),
        actionButton(paste0(p, "_padd"), "Add this pair", icon = icon("plus"), class = "btn-light btn-sm w-100")),
      uiOutput(paste0(p, "_pairs_ui"))),
    tags$label(class = "control-label", "Test", tip_options(TEST_CH, STAT_TIPS$test)),
    selectInput(paste0(p, "_test"), NULL, TEST_CH), opt_hints(p, "test", STAT_TIPS$test),
    tags$label(class = "control-label", "Multiple-testing correction", tip_options(ADJ_CH, STAT_TIPS$adj)),
    selectInput(paste0(p, "_adj"), NULL, ADJ_CH), opt_hints(p, "adj", STAT_TIPS$adj),
    selectInput(paste0(p, "_lab"), lab("Show as", "p-values: the adjusted p-value written out. Stars: * p < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001, ns not significant."), c("p-values" = "p", "Stars" = "stars")),
    input_switch(paste0(p, "_hidens"), lab("Hide non-significant", "Only results with adjusted p < 0.05 are drawn. The table keeps all of them."), hide),
    num_step(paste0(p, "_statsize"), "Label size", 3.5, 2, 8, 0.25))
}

style_panel <- function(p, h = 480, legend = TRUE, ...) accordion_panel("Style", value = "style", icon = icon("brush"),
  div(class = "mcx-grid2",
    num_step(paste0(p, "_pt"), "Point size", 1.2, 0.2, 5, 0.2),
    num_step(paste0(p, "_alpha"), "Point opacity", 0.6, 0.1, 1, 0.1),
    num_step(paste0(p, "_txt"), "Text size", 12, 8, 24, 1),
    num_step(paste0(p, "_height"), "Height (px)", h, 260, 1800, 40)),
  selectInput(paste0(p, "_angle"), "Axis labels", c("Automatic" = "auto", "Flat" = "0", "Tilted" = "45", "Vertical" = "90")),
  input_switch(paste0(p, "_legend"), "Legend", legend),
  input_switch(paste0(p, "_mean"), lab("Line at the overall mean", "Box, violin, bar and strip plots."), FALSE),
  input_switch(paste0(p, "_values"), lab("Numbers on heatmap tiles", "Heatmap only."), FALSE), ...)
more <- function(..., id = NULL) accordion(id = id, open = FALSE, multiple = TRUE, class = "mcx-more", ...)

# "How to use": step cards, an example drawn as a flow of chips, and longer notes folded away
HELP_VISUAL <- list(annotate = function() {
  cols <- c("#B22222", "#E8A2A2", "#F5F5F5", "#9DB0D8", "#24408E")
  cell <- function(v) span(class = "c", style = paste0("background:", cols[v]))
  set.seed(3); m <- matrix(sample(2:5, 6 * 12, TRUE), 6); m[1:3, 5:7] <- 1
  div(class = "mcx-hvis",
    div(class = "mcx-hvis-hm", div(class = "bar", lapply(1:12, function(j) span(class = "c", style = paste0("background:", if (j %in% 5:7) "#8C6D62" else "#E5C35A")))),
        div(class = "mcx-g", lapply(1:6, function(i) div(class = "r", lapply(1:12, function(j) cell(m[i, j]))))), div(class = "sel"), div(class = "cap", "Click B22")),
    div(class = "mcx-hvis-arrow", icon("arrow-right")),
    div(class = "mcx-hvis-panel", tags$b("B22"), div("84% Mono/Mac"), div("M12 C1QA C1QB"), div(class = "nm", icon("tag"), " Resident macrophages")))
})
howto <- function(key) {
  h <- HELP[[key]]; if (is.null(h)) return(NULL); id <- paste0("help_", key)
  chip <- function(x) { p <- strsplit(x, "|", fixed = TRUE)[[1]]; span(class = "mcx-hchip", span(class = "k", p[1]), span(class = "v", HTML(p[2]))) }
  tagList(
    tags$button(type = "button", class = "btn mcx-howto", `data-bs-toggle` = "modal", `data-bs-target` = paste0("#", id), icon("circle-question"), "How to use"),
    div(class = "modal fade mcx-howto-modal", id = id, tabindex = "-1", `aria-hidden` = "true",
      div(class = "modal-dialog modal-lg modal-dialog-scrollable", div(class = "modal-content",
        div(class = "modal-header", div(h5(class = "modal-title", h$title), div(class = "mcx-howto-what", h$what)),
            tags$button(type = "button", class = "btn-close", `data-bs-dismiss` = "modal", `aria-label` = "Close")),
        div(class = "modal-body",
          div(class = "mcx-hsteps", lapply(seq_along(h$steps), function(i) { st <- h$steps[[i]]
            div(class = "mcx-hstep", div(class = "top", span(class = "n", i), icon(st[1])), div(class = "t", st[2]), div(class = "x", HTML(st[3]))) })),
          if (!is.null(h$visual)) HELP_VISUAL[[h$visual]](),
          if (!is.null(h$example)) div(class = "mcx-hexample", span(class = "lbl", "Example"),
            div(class = "flow", lapply(seq_along(h$example$flow), function(i) tagList(if (i > 1) span(class = "arr", icon("chevron-right")), chip(h$example$flow[i])))),
            if (!is.null(h$example$result)) div(class = "res", icon("circle-check"), HTML(h$example$result))),
          if (!is.null(h$manual)) div(class = "mcx-manual", lapply(h$manual, function(sec)
            div(class = "mcx-msec", div(class = "mcx-msec-h", icon(sec[[1]]), sec[[2]]),
                lapply(sec[[3]], function(it) if (!nzchar(it[1])) div(class = "mcx-mrow plain", HTML(it[2]))
                       else div(class = "mcx-mrow", span(class = "mcx-mctl", it[1]), span(class = "mcx-mtxt", HTML(it[2])))))))))))
  )
}
page_head <- function(title, sub, start = NULL, help = NULL) div(class = "mcx-head",
  div(class = "mcx-head-text", h2(title), p(sub), if (!is.null(start)) p(class = "mcx-start", tags$b("Start here: "), start)),
  if (!is.null(help)) howto(help))
empty_state <- function(ic, title, text, goto = NULL) div(class = "mcx-empty", div(class = "i", icon(ic)), h5(title), p(text),
  if (!is.null(goto)) tags$a(href = "#", class = "btn btn-primary", onclick = sprintf("Shiny.setInputValue('goto','%s',{priority:'event'});return false;", goto), "Open a dataset"))
needs_data <- function(...) tagList(
  conditionalPanel("!output.has_data", empty_state("folder-open", "No dataset open", "A dataset is opened on the Home page.", goto = "home")),
  conditionalPanel("output.has_data", ...))
sidebar_page <- function(title, sub, side, ..., start = NULL, help = NULL, width = 280, top = NULL)
  div(class = "mcx-page", page_head(title, sub, start, help),
      needs_data(top, layout_sidebar(sidebar = sidebar(width = width, open = "desktop", side), border = TRUE, fillable = FALSE, ...)))

dl_bar <- function(id, w = 8, h = 6) {
  ns <- NS(id)
  div(class = "dl-bar",
      span(class = "dl-title", icon("download"), "Figure"),
      selectInput(ns("fmt"), NULL, c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg"), width = "74px", selectize = FALSE),
      numericInput(ns("w"), NULL, w, min = 1, max = 80, step = 0.5, width = "60px"), span("x"),
      numericInput(ns("h"), NULL, h, min = 1, max = 80, step = 0.5, width = "60px"), span("inches"),
      downloadButton(ns("go"), "Download", class = "btn-sm btn-outline-primary"))
}
dl_server <- function(id, draw, fname) moduleServer(id, function(input, output, session) {
  output$go <- downloadHandler(filename = function() paste0(fname(), ".", input$fmt),
                               content = function(file) save_plot(file, draw(), input$w, input$h, input$fmt))
})
plot_card <- function(id, dl_id, w = 8, h = 6, header = NULL, ..., click = NULL, brush = NULL, below = NULL)
  card(full_screen = TRUE, class = "mcx-plotcard", if (!is.null(header)) card_header(header),
       card_body(..., plotOutput(id, height = "auto", click = click, brush = brush), below), card_footer(dl_bar(dl_id, w, h)))
# table header with a hover explanation per column
dt_head <- function(cols, tips) htmltools::withTags(table(class = "display", thead(tr(lapply(cols, function(n)
  th(title = unname(tips[n]) %|na|% NULL, if (!is.na(tips[n])) span(n, span(class = "mcx-tip", "?")) else n))))))
# tables live in collapsible cards so long pages stay short
table_card <- function(title, id, ..., help = NULL, tools = NULL, open = TRUE, icon_name = "table")
  fold(span(title, if (!is.null(help)) tip(help)), ..., DTOutput(id), icon_name = icon_name, open = open, tools = tools)
