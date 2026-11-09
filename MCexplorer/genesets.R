# =============================================================================
#  MCexplorer  |  genesets.R
#  Published gene sets for signatures: MSigDB collections through the msigdbr package
#  (optional; install.packages("msigdbr"), and for msigdbr >= 10 also the msigdbdf data package),
#  or any .gmt file (the format MSigDB, Enrichr and GSEA download in).
# =============================================================================

GENESET_COLLECTIONS <- c(
  "Hallmark (50 well-defined processes)" = "H",
  "GO biological process"                = "GO:BP",
  "Reactome pathways"                    = "CP:REACTOME",
  "KEGG pathways (legacy)"               = "CP:KEGG_LEGACY",
  "WikiPathways"                         = "CP:WIKIPATHWAYS",
  "Cell type signatures (C8)"            = "C8")
GENESET_SPEC <- list(H = c("H", NA), "GO:BP" = c("C5", "GO:BP"), "CP:REACTOME" = c("C2", "CP:REACTOME"),
                     "CP:KEGG_LEGACY" = c("C2", "CP:KEGG_LEGACY"), "CP:WIKIPATHWAYS" = c("C2", "CP:WIKIPATHWAYS"), C8 = c("C8", NA))

genesets_ok <- function() requireNamespace("msigdbr", quietly = TRUE)

# human gene symbols are upper case (CD3D), mouse ones are not (Cd3d)
guess_species <- function(genes) {
  g <- head(genes[grepl("[A-Za-z]", genes)], 2000)
  if (mean(g == toupper(g)) > 0.6) "Homo sapiens" else "Mus musculus"
}

# one MSigDB collection as a named list of gene vectors; works with old (category) and new (collection) msigdbr
load_genesets <- function(code, species) {
  if (!genesets_ok()) stop("The gene-set library needs the msigdbr package: install.packages(\"msigdbr\")")
  sp <- GENESET_SPEC[[code]]; f <- msigdbr::msigdbr; new_api <- "collection" %in% names(formals(f))
  subs <- if (identical(code, "CP:KEGG_LEGACY")) c("CP:KEGG_LEGACY", "CP:KEGG") else sp[2]
  df <- NULL; err <- NULL
  for (s in subs) {
    sub <- if (is.na(s)) NULL else s
    df <- tryCatch(if (new_api) f(species = species, collection = sp[1], subcollection = sub) else f(species = species, category = sp[1], subcategory = sub),
                   error = function(e) { err <<- conditionMessage(e); NULL })
    if (!is.null(df) && nrow(df)) break
  }
  if (is.null(df) || !nrow(df)) stop(if (!is.null(err) && grepl("msigdbdf", err)) "msigdbr 10 or newer keeps the gene sets in a second package: install.packages(\"msigdbdf\", repos = \"https://igordot.r-universe.dev\")"
                                     else paste("This collection could not be loaded.", err %||% ""))
  gs <- split(as.character(df$gene_symbol), df$gs_name)
  lapply(gs, unique)
}

read_gmt <- function(path) {
  l <- strsplit(readLines(path, warn = FALSE), "\t")
  l <- l[lengths(l) >= 3]
  if (!length(l)) stop("No gene sets found. A .gmt file has one set per line: name, description, then genes, separated by tabs.")
  setNames(lapply(l, function(x) { g <- trimws(x[-(1:2)]); unique(g[nzchar(g)]) }), vapply(l, `[`, "", 1))
}

# "HALLMARK_INTERFERON_GAMMA_RESPONSE" -> "Interferon gamma response (Hallmark)"
nice_set_name <- function(x) {
  pre <- sub("_.*$", "", x); rest <- sub("^[^_]+_", "", x)
  src <- c(HALLMARK = "Hallmark", GOBP = "GO", REACTOME = "Reactome", KEGG = "KEGG", WP = "WikiPathways")[pre]
  rest <- gsub("_", " ", tolower(rest)); rest <- paste0(toupper(substr(rest, 1, 1)), substr(rest, 2, nchar(rest)))
  ifelse(is.na(src), x, paste0(rest, " (", src, ")"))
}
