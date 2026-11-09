# =============================================================================
#  MCexplorer  |  server.R
# =============================================================================

RESTORE <- c(grp_by = "select", grp_k = "numeric", grp_seed = "numeric", hm_mods = "multi", hm_view = "radio", hm_bars = "multi", hm_height = "slider", hm_modsel = "radio",
  hm_sep = "switch", hm_glab = "switch", hm_lim = "slider", hm_low = "colour", hm_mid = "colour", hm_high = "colour",
  g_src = "select", g_top_n = "numeric", g_var_n = "numeric", g_custom = "text", h2_view = "radio", h2_bars = "multi", h2_modsel = "radio", h2_lim = "slider", h2_gpal = "select",
  h2_genez = "switch", h2_sep = "switch", h2_gprop = "slider", h2_height = "slider",
  pl_x = "select", pl_what = "radio", pl_metric = "select", pl_split = "select", pl_type = "radio",
  pl_points = "switch", pl_stats = "switch", pl_scale = "switch", pl_cmp = "select", pl_test = "select", pl_adj = "select", pl_lab = "select", pl_hidens = "switch",
  pl_pt = "slider", pl_alpha = "slider", pl_txt = "slider", pl_angle = "select", pl_legend = "switch", pl_mean = "switch", pl_height = "slider",
  cd_mode = "radio", cd_cond = "select", cd_unit = "select", cd_grp_by = "select", cd_stack_by = "radio", cd_cells = "multi", cd_what = "radio",
  cd_type = "radio", cd_points = "switch", cd_stats = "switch", cd_scale = "switch", cd_cmp = "select", cd_test = "select", cd_adj = "select", cd_lab = "select",
  cd_hidens = "switch", cd_min = "numeric", cd_pt = "slider", cd_alpha = "slider", cd_txt = "slider", cd_ncol = "numeric", cd_height = "slider",
  um_by = "select", um_labels = "switch", um_pt = "slider", um_height = "slider", um_pal = "select")

server <- function(input, output, session) {

  D <- reactiveVal(NULL)       # core dataset
  CELLS <- reactiveVal(NULL)   # cell information summarised per metacell
  UMAP <- reactiveVal(NULL)    # list(xy, source)
  rv <- reactiveValues(scan = NULL, up = NULL, ann = NULL, hist = list(), log = NULL, colors = character(0),
                       job = NULL, cell_status = list(state = "none"), vg = NULL, session_msg = NULL, ann_msg = NULL,
                       umap_msg = NULL, gops = list(), split_sig = "", sel = character(0), tbl_sync = NULL,
                       ce_levels = NULL, sigs = list(), sig_target = NULL, sig_select = NULL, sig_edit = NULL,
                       labels_applied = FALSE, pairs = list(pl = character(0), cd = character(0)), pending = list(),
                       figs_seen = character(0), tables_seen = character(0), open_msg = NULL, wizard = FALSE)
  colmem <- new.env(); cache_env <- new.env(); gs_env <- new.env(); score_env <- new.env()
  figs <- new.env(); tables <- new.env()

  disable <- function(ids, state = TRUE, reason = NULL) session$sendCustomMessage("mcx-disable", list(ids = as.list(ids), disable = state, reason = reason))
  disable_choice <- function(id, value, state = TRUE, reason = NULL) session$sendCustomMessage("mcx-disable-choice", list(id = id, value = value, disable = state, reason = reason))
  go <- function(tab) nav_select("nav", tab)
  observeEvent(input$goto, go(input$goto))
  note <- function(...) showNotification(paste0(...), type = "message", duration = 4)
  # every table goes through plain columns (no arrays, lists or factors), which DataTables can always send
  datatable <- function(data, ...) {
    data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
    for (j in seq_along(data)) { x <- data[[j]]
      data[[j]] <- if (is.list(x)) vapply(x, function(z) paste(unlist(z), collapse = ", "), "") else if (is.factor(x)) as.character(x)
                   else if (is.numeric(x) || is.logical(x)) as.vector(unname(x)) else as.character(unname(x)) }
    DT::datatable(data, ...)
  }
  # numeric input with a fallback while the box is empty or being typed in
  iv <- function(id, d) { v <- suppressWarnings(as.numeric(input[[id]])); if (!length(v) || is.na(v)) d else v }
  # a warning says what went wrong and, as its last argument, what to do next
  warn <- function(..., todo = NULL) {
    showNotification(tagList(paste0(..., collapse = ""), if (!is.null(todo)) span(class = "mcx-do", icon("arrow-right"), " ", todo)), type = "warning", duration = 7)
  }

  output$has_data  <- reactive(!is.null(D()))
  output$has_cells <- reactive(length(CELLS()$comp) > 0 && isTRUE(CELLS()$cell_level))
  output$scan_multi <- reactive(length(rv$scan$sets) > 1)
  output$show_wizard <- reactive(is.null(D()) || isTRUE(rv$wizard))
  observeEvent(input$show_open, { rv$wizard <- !isTRUE(rv$wizard) })
  output$scan_has_cells <- reactive(length(rv$scan$cell_cand) > 0)
  output$scan_has_umap <- reactive(length(rv$scan$umap_cand) > 0)
  for (o in c("has_data", "has_cells", "show_wizard", "scan_multi", "scan_has_cells", "scan_has_umap")) outputOptions(output, o, suspendWhenHidden = FALSE)

  reg_fig <- function(key, label, draw, w, h) {
    assign(key, list(draw = draw, w = w, h = h, label = label), envir = figs)
    if (!key %in% isolate(rv$figs_seen)) rv$figs_seen <- c(isolate(rv$figs_seen), key)
  }
  reg_table <- function(key, x) {
    assign(key, x, envir = tables)
    if (!key %in% isolate(rv$tables_seen)) rv$tables_seen <- c(isolate(rv$tables_seen), key)
  }

  # ============================================================ colours
  pal_for <- function(lv) {
    lv <- as.character(lv); user <- rv$colors
    vapply(lv, function(x) {
      if (!is.na(user[x])) return(unname(user[x]))
      if (is.na(x) || x %in% NO_LABEL || x == "Not selected") return(GREY)
      if (x == "Selected") return(ACCENT)
      if (!exists(x, envir = colmem, inherits = FALSE)) {
        if (grepl("^[BN][0-9]+(\\.[0-9]+)?$", x)) assign(x, spread_col(as.numeric(sub("^[BN]", "", x)) * 10), envir = colmem)   # blocks / clusters
        else { n <- sum(!grepl("^[BN][0-9]+(\\.[0-9]+)?$", ls(colmem, all.names = TRUE)))
               assign(x, if (n < length(PAL)) PAL[n + 1] else spread_col(n), envir = colmem) }
      }
      get(x, envir = colmem)
    }, "", USE.NAMES = TRUE)
  }
  heat_cols <- reactive(c(col_or(input$hm_low, "#000080"), col_or(input$hm_mid, "#FFFFFF"), col_or(input$hm_high, "#B22222")))
  observeEvent(input$mcx_colour, { x <- input$mcx_colour; if (valid_col(x$value %||% "")) colourpicker::updateColourInput(session, x$id, value = toupper(x$value)) })
  # colours window: an optional heatmap colour scale (scale = list(ids, defaults, extra UI)) and the group colours
  open_colors <- function(levels, title = "Colours", scale = NULL) {
    levels <- unique(as.character(levels[!is.na(levels)]))
    if (!length(levels) && is.null(scale)) return(warn("Nothing to colour here yet.", todo = "Colours can be set once groups are shown in the plot."))
    if (length(levels) > 60) { note("Showing the first 60 groups."); levels <- head(levels, 60) }
    rv$ce_levels <- levels; pal <- pal_for(levels)
    scale_ui <- if (!is.null(scale)) div(class = "mcx-colscale",
      tags$label(class = "control-label", lab("Heatmap colour scale", "Low, middle and high end of the module score colours.")),
      div(class = "mcx-colpair", lapply(1:3, function(i) colour_field(scale$ids[i], c("Low", "Middle", "High")[i], col_or(input[[scale$ids[i]]], scale$def[i])))),
      scale$extra)
    showModal(modalDialog(title = title, size = "m", easyClose = TRUE, scale_ui,
      if (length(levels)) tagList(if (!is.null(scale)) tags$label(class = "control-label mt-2", "Label and group colours"),
      p(class = "mcx-hint", "Each colour can be typed as a hex code (#1F77B4) or a name (steelblue), picked from the square, the colour wheel ",
        icon("circle-half-stroke"), ", or anywhere on screen with the eyedropper ", icon("eye-dropper"), "."),
      div(class = "mcx-inline mb-2",
        div(class = "flex-grow-1", selectInput("ce_pal", "Fill all from a palette", names(PALETTES), width = "100%", selectize = FALSE)),
        actionButton("ce_pal_go", "Fill", class = "btn-light")),
      div(class = "mcx-colgrid", lapply(seq_along(levels), function(i) div(class = "mcx-colrow", span(class = "mcx-colname", levels[i]),
        colour_field(paste0("ce_", i), NULL, hex6(pal[[i]]))))),
      textAreaInput("ce_paste", lab("Paste several at once", "One per line: name = #hexcode"), rows = 2, width = "100%",
                    placeholder = "CD8 T = #4E79A7\nMono/Mac = #E15759")),
      footer = tagList(actionButton("ce_reset", "Back to defaults", class = "btn-light"), modalButton("Cancel"),
                       actionButton("ce_apply", "Apply", class = "btn-primary"))))
  }
  observeEvent(input$ce_pal_go, { lv <- rv$ce_levels; cols <- palette_cols(input$ce_pal, length(lv))
    for (i in seq_along(lv)) colourpicker::updateColourInput(session, paste0("ce_", i), value = hex6(cols[i])) })
  observeEvent(input$ce_apply, {
    lv <- rv$ce_levels; new <- rv$colors
    for (i in seq_along(lv)) { v <- trimws(input[[paste0("ce_", i)]] %||% ""); if (valid_col(v)) new[lv[i]] <- hex6(v) }
    for (ln in strsplit(input$ce_paste %||% "", "\n")[[1]]) {
      kv <- trimws(strsplit(ln, "=|\t")[[1]]); if (length(kv) == 2 && nzchar(kv[1]) && valid_col(kv[2])) new[kv[1]] <- hex6(kv[2])
    }
    rv$colors <- new; removeModal()
  })
  observeEvent(input$ce_reset, { rv$colors <- rv$colors[setdiff(names(rv$colors), rv$ce_levels)]; removeModal() })

  # ============================================================ values per metacell
  meta_choices <- function() {
    cl <- CELLS(); l <- names(cl$labels); m <- names(cl$comp)
    c(if (length(l)) setNames(paste0("lab:", l), pretty_name(l)), if (length(m)) setNames(paste0("maj:", m), paste("Main", pretty_name(m))))
  }
  primary_label <- reactive({
    cl <- CELLS(); if (is.null(cl)) return(NULL)
    bad <- "sample|donor|patient|orig|batch|condition|disease|site|tissue|phase|block|module|metacell|treatment|status"
    cand <- c(if (length(cl$labels)) setNames(paste0("lab:", names(cl$labels)), names(cl$labels)),
              if (length(cl$comp)) setNames(paste0("maj:", names(cl$comp)), names(cl$comp)))
    if (!length(cand)) return(NULL)
    nm <- names(cand); ok <- grepl("label|type|annot|lineage|ident|class", nm, ignore.case = TRUE) & !grepl(bad, nm, ignore.case = TRUE)
    if (any(ok)) { k <- cand[ok]; return(unname(c(k[grepl("broad", names(k), ignore.case = TRUE)], k)[1])) }
    k <- cand[startsWith(cand, "lab:") & !grepl(bad, nm, ignore.case = TRUE)]
    if (length(k)) unname(k[1]) else NULL
  })
  comp_cols <- reactive(names(CELLS()$comp) %||% character(0))
  guess_cond <- reactive({
    cc <- comp_cols(); if (!length(cc)) return(NULL)
    g <- grep("condition|disease|status|group|treatment|arm", cc, value = TRUE, ignore.case = TRUE)
    if (length(g)) return(g[1])
    nl <- vapply(cc, function(c) ncol(CELLS()$comp[[c]]), 1); (cc[nl >= 2 & nl <= 8] %||% cc)[1]
  })
  unit_candidates <- function(cond) {
    cl <- CELLS(); cc <- setdiff(comp_cols(), cond); if (!length(cc)) return(character(0))
    pur <- vapply(cc, function(u) { t <- cl$pairs[[u]][[cond]]; if (is.null(t)) return(0); t <- t[rowSums(t) > 0, , drop = FALSE]; min(apply(t, 1, max) / rowSums(t)) }, 1)
    nl <- vapply(cc, function(u) ncol(cl$comp[[u]]), 1)
    cc[order(-((pur >= 0.95) * 1000 + grepl("sample|donor|patient|orig", cc, ignore.case = TRUE) * 100 + nl))]
  }

  mv <- function(key, mcs = NULL) {
    b <- D(); cl <- CELLS(); if (is.null(mcs)) mcs <- b$mcs
    if (key == "ann") return(rv$ann$annotation[match(mcs, rv$ann$metacell)])
    if (key == "block") return(paste0("B", b$block[mcs]))
    if (key == "grp") return(unname(grouping()$id[mcs]))
    if (key == "sel") return(ifelse(mcs %in% rv$sel, "Selected", "Not selected"))
    type <- sub(":.*$", "", key); col <- sub("^[^:]*:", "", key); i <- match(mcs, b$mcs)
    switch(type,
      lab = cl$labels[[col]][i], maj = cl$comp_major[[col]][i], pur = cl$comp_purity[[col]][i],
      num = switch(col, n_genes = b$nfeat[i], n_cells = cl$n_cells[i], total_umis = cl$total_umis[i], cl$nums[[col]][i]),
      NULL)
  }
  key_label <- function(key) {
    switch(key, ann = return("Names given"), block = return("Pipeline block"), grp = return("Annotate block"), sel = return("Selection"))
    col <- sub("^[^:]*:", "", key)
    switch(sub(":.*$", "", key), lab = pretty_name(col), maj = paste("Main", pretty_name(col)), pur = paste("Share of main", pretty_name(col)),
           num = switch(col, n_genes = "Genes detected", n_cells = "Cells", total_umis = "UMIs", paste("Mean", pretty_name(col))), col)
  }
  level_order <- function(key, v) {
    v <- unique(as.character(v[!is.na(v)]))
    if (key == "block") return(intersect(paste0("B", seq_len(D()$nb)), v))
    if (key == "grp") return(intersect(grouping()$levels, v))
    if (key == "sel") return(intersect(c("Selected", "Not selected"), v))
    if (key == "ann") return(c(sort(setdiff(v, NO_LABEL)), intersect(v, NO_LABEL)))
    col <- sub("^[^:]*:", "", key); lv <- CELLS()$levels[[col]]
    if (!is.null(lv)) return(c(intersect(lv, v), setdiff(v, lv)))
    mixed_sort(v)
  }
  group_choices <- function(sel = TRUE, meta = TRUE) {
    out <- list("Annotate page" = c("Names given" = "ann", "Annotate blocks" = "grp", if (sel) c("Current selection" = "sel")),
                "Pipeline" = c("Pipeline blocks" = "block"))
    mc <- if (meta) meta_choices() else NULL
    if (length(mc)) out[["Cell metadata"]] <- mc
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
  group_major <- function(key, f, pct = FALSE) {
    v <- mv(key); if (is.null(v)) return(rep(NA_character_, nlevels(f)))
    v <- as.character(v); keep <- !is.na(v) & !is.na(f)
    t <- table(f[keep], v[keep]); if (!ncol(t)) return(rep(NA_character_, nlevels(f)))
    top <- colnames(t)[max.col(t, "first")]; mx <- apply(t, 1, max); rs <- rowSums(t)
    out <- ifelse(rs == 0, NA_character_, if (pct) ifelse(mx < rs, paste0(top, " (", round(100 * mx / pmax(rs, 1)), "%)"), top) else top)
    unname(out[levels(f)])
  }
  # main value of a label in each group, as a share of the group's cells (not of its metacells)
  group_cell_major <- function(key, f) {
    cl <- CELLS(); type <- sub(":.*$", "", key); col <- sub("^[^:]*:", "", key)
    if (type == "maj" && !is.null(cl$comp[[col]])) tab <- rowsum(cl$comp[[col]], f, reorder = TRUE)
    else { v <- as.character(mv(key)); if (!length(v)) return(rep(NA_character_, nlevels(f))); v[is.na(v)] <- "NA"
      w <- if (!is.null(cl$n_cells)) cl$n_cells else rep(1, length(v)); tab <- unclass(tapply(w, list(f, v), sum)); tab[is.na(tab)] <- 0 }
    tab <- tab[match(levels(f), rownames(tab)), , drop = FALSE]
    top <- colnames(tab)[max.col(tab, "first")]; sh <- apply(tab, 1, max) / pmax(rowSums(tab), 1)
    ifelse(rowSums(tab) == 0, NA_character_, paste0(top, " (", round(100 * sh), "% of cells)"))
  }
  # "B22" -> "B22 · Mono/Mac": blocks and clusters carry their main label so they can be read
  disp_labels <- function(key, lv) {
    out <- setNames(as.character(lv), lv)
    by <- if (key == "grp") grouping()$by else key
    if (!key %in% c("block", "grp") || !by %in% c("block", "cluster") || length(lv) > 400) return(out)
    src <- primary_label() %||% (if (any(!rv$ann$annotation %in% NO_LABEL)) "ann")
    if (is.null(src)) return(out)
    f <- factor(mv(key), levels = lv); m <- group_major(src, f)
    out[] <- ifelse(is.na(m), lv, paste0(lv, " · ", m))
    out
  }
  mod_choice_df <- function() {
    b <- D(); o <- order(b$mods$num)
    data.frame(value = b$mods$id[o], label = b$mods$lab[o], short = b$mods$short[o], genes = b$mods$genes[o],
               first = vapply(b$mod_genes[b$mods$id[o]], function(g) paste(head(g, 4), collapse = " "), ""), stringsAsFactors = FALSE)
  }
  # chosen modules show as "M12"; the drop-down also shows each module's first genes
  MOD_RENDER <- function(short) I(paste0("{item: function(d, esc) { return '<div class=\"item\">' + esc(d.", if (short) "short" else "label", ") + '</div>'; },
    option: function(d, esc) { return '<div class=\"option\"><b>' + esc(d.short) + '</b><span class=\"mcx-opt-sub\">' + esc(d.first) + '</span></div>'; }}"))
  set_mod_choices <- function(id, selected = NULL) {
    df <- mod_choice_df()
    updateSelectizeInput(session, id, choices = df, selected = selected %||% df$value[1], server = TRUE,
                         options = list(searchField = c("label", "genes"), placeholder = "Search a module number or a gene",
                                        render = MOD_RENDER(id %in% c("pl_mod", "cd_mod", "h2_mods", "hm_mods"))))
  }
  set_gene_choices <- function(id, selected) updateSelectizeInput(session, id, choices = rownames(D()$counts), selected = selected, server = TRUE)
  # a gene list as plain text (not buttons), with an optional link that copies the whole list
  genes_ui <- function(g, n = 60, copy = FALSE, label = "Genes:") {
    div(class = "mcx-genetext", span(class = "lbl", label), paste(head(g, n), collapse = ", "),
        if (length(g) > n) span(class = "text-muted", paste0(" and ", length(g) - n, " more")),
        if (copy) tags$button(type = "button", class = "mcx-copy", `data-copy` = paste(g, collapse = ", "), icon("copy"), " Copy all"))
  }

  # ------------------------------------------------------------ draggable lists
  rank_set <- function(id, values, labels = values, on = TRUE, colors = NULL) {
    on <- rep_len(on, length(values))
    items <- lapply(seq_along(values), function(i) list(value = values[i], label = unname(labels[i]) %||% values[i], on = isTRUE(on[i]),
                                                        color = if (is.null(colors)) "" else unname(colors[i])))
    session$sendInputMessage(id, list(items = items))
  }
  group_levels <- function(key) { v <- mv(key); if (is.null(v)) character(0) else level_order(key, v) }
  rank_prev <- new.env()
  fill_groups_rank <- function(id, key, ticked = NULL, default_n = Inf, max_items = 500) {
    lv <- head(group_levels(key), max_items); dl <- disp_labels(key, lv)
    pend <- isolate(rv$pending[[id]]); if (!is.null(pend)) { ticked <- pend; rv$pending[[id]] <- NULL }
    prev <- rank_prev[[id]]; cur <- isolate(input[[id]])
    if (is.null(ticked) && identical(prev$key, key) && length(intersect(cur, lv))) ticked <- c(intersect(cur, lv), setdiff(lv, prev$lv))
    if (is.null(ticked)) { real <- setdiff(lv, NO_LABEL); ticked <- head(if (length(real)) real else lv, default_n) }
    assign(id, list(key = key, lv = lv), envir = rank_prev)
    lv <- c(intersect(ticked, lv), setdiff(lv, ticked))
    rank_set(id, lv, dl[lv], lv %in% ticked, pal_for(lv))
  }
  # select a grouping key in `key_input` and tick `ticked` in rank list `id`, whether or not the key changes
  set_groups <- function(key_input, key, id, ticked) {
    if (identical(isolate(input[[key_input]]), key)) fill_groups_rank(id, key, ticked)
    else { rv$pending[[id]] <- ticked; updateSelectInput(session, key_input, selected = key) }
  }

  # ============================================================ HOME: opening data
  do_scan <- function(dir) {
    dir <- sub("/+$", "", trimws(dir)); if (!nzchar(dir)) return()
    if (!dir.exists(dir)) { rv$scan <- list(error = paste0("Folder not found: ", dir)); return() }
    s <- scan_folder(dir); rv$scan <- s
    if (length(s$sets)) updateSelectInput(session, "set", choices = names(s$sets), selected = names(s$sets)[1])
    if (length(s$cell_cand)) updateSelectInput(session, "cell_file", choices = setNames(s$cell_cand, basename(s$cell_cand)),
                                               selected = guess_cell_file(s$cell_cand, names(s$sets)[1] %||% ""))
    if (length(s$umap_cand)) updateSelectInput(session, "umap_file", choices = setNames(s$umap_cand, basename(s$umap_cand)))
    if (!length(s$cell_cand)) updateRadioButtons(session, "cell_mode", selected = if (nzchar(input$cell_path %||% "")) "file" else "none")
    if (length(s$umap_cand)) updateRadioButtons(session, "umap_mode", selected = "file")
    if (length(s$sets) == 1 && identical(dir, sub("/+$", "", DATA_PATH)) && is.null(D())) open_folder()
  }
  folder_d <- debounce(reactive(trimws(input$folder %||% "")), 700)   # files are found as soon as the path is typed
  observeEvent(folder_d(), { if (nzchar(folder_d())) do_scan(folder_d()) }, ignoreInit = TRUE)
  observeEvent(TRUE, { if (nzchar(DATA_PATH)) do_scan(DATA_PATH) }, once = TRUE)
  observeEvent(input$set, { s <- rv$scan; req(length(s$cell_cand)); updateSelectInput(session, "cell_file", selected = guess_cell_file(s$cell_cand, input$set)) }, ignoreInit = TRUE)

  file_li <- function(ok, label, path = NULL, optional = FALSE)
    tags$li(span(class = if (ok) "ok" else if (optional) "opt" else "no", icon(if (ok) "circle-check" else if (optional) "circle-minus" else "circle-xmark")),
            span(strong(label), if (!is.null(path) && !is.na(path)) tagList(" ", code(basename(path)))))
  output$scan_msg <- renderUI({
    s <- rv$scan; if (is.null(s)) return(NULL)
    if (!is.null(s$error)) return(div(class = "alert alert-warning py-2", icon("triangle-exclamation"), " ", s$error))
    if (!length(s$sets)) {
      miss <- if (length(s$incomplete)) paste0(" Missing: ", paste(ROLE_NAMES[setdiff(names(ROLE_KEYS), names(s$incomplete[[1]]))], collapse = ", "), ".") else ""
      return(div(class = "alert alert-warning py-2", icon("triangle-exclamation"), " The four needed .rd files were not all found in this folder.", miss))
    }
    NULL
  })
  output$scan_files <- renderUI({
    s <- rv$scan; req(length(s$sets)); f <- s$sets[[input$set %||% names(s$sets)[1]]]; req(f)
    tags$ul(class = "mcx-files", lapply(names(ROLE_KEYS), function(r) file_li(!is.na(f[r]), ROLE_NAMES[[r]], f[r])))
  })
  observeEvent(input$up_req, {
    fi <- input$up_req; out <- rv$up %||% list(objs = list(), names = list(), unknown = character(0)); out$unknown <- character(0)
    withProgress(message = "Reading uploaded files", value = 0, {
      for (i in seq_len(nrow(fi))) {
        nm <- fi$name[i]; incProgress(1 / nrow(fi), detail = nm)
        e <- new.env(); objs <- tryCatch(load(fi$datapath[i], envir = e), error = function(err) character(0)); found <- FALSE
        for (o in objs) { r <- rd_role(get(o, envir = e)); if (!is.na(r)) { out$objs[[r]] <- get(o, envir = e); out$names[[r]] <- nm; found <- TRUE } }
        if (!found) out$unknown <- c(out$unknown, nm)
      }
    })
    rv$up <- out
  })
  output$up_files_ui <- renderUI({
    u <- rv$up; if (is.null(u)) return(NULL)
    tagList(tags$ul(class = "mcx-files", lapply(names(ROLE_KEYS), function(r) file_li(r %in% names(u$objs), ROLE_NAMES[[r]], u$names[[r]]))),
            if (length(u$unknown)) div(class = "alert alert-warning py-2", "Not recognised: ", paste(u$unknown, collapse = ", ")))
  })
  open_folder <- function() {
    s <- rv$scan
    if (is.null(s) || !length(s$sets)) return(warn("No module result files found yet.", todo = "Paste the folder path in step 2 and click Find files."))
    stem <- input$set %||% names(s$sets)[1]; if (!stem %in% names(s$sets)) stem <- names(s$sets)[1]
    cell <- if (identical(input$cell_mode, "file")) { p <- trimws(input$cell_path %||% ""); if (nzchar(p)) p else (input$cell_file %||% guess_cell_file(s$cell_cand, stem) %||% "") } else ""
    umap <- if (identical(input$umap_mode, "file")) { p <- trimws(input$umap_path %||% ""); if (nzchar(p)) p else (input$umap_file %||% "") } else ""
    if (nzchar(cell) && !file.exists(cell)) return(warn("Cell information file not found: ", cell, todo = "Check the path in step 3, or choose Skip."))
    if (nzchar(umap) && !file.exists(umap)) return(warn("UMAP file not found: ", umap, todo = "Check the path in step 4, or choose Compute from the module scores."))
    open_dataset(as.list(s$sets[[stem]]), stem, cache_path(s$dir, stem), cell, umap)
  }
  observeEvent(input$open_go, {
    if (identical(input$src, "upload")) {
      u <- rv$up
      if (is.null(u) || !all(names(ROLE_KEYS) %in% names(u$objs))) return(warn("The four module result files are not all uploaded.", todo = "Choose all four .rd files in step 2."))
      stem <- gsub("_+", "_", gsub("(genelist|metadata|counts|hmlist)", "", tolower(sub("\\.[^.]+$", "", u$names$counts))))
      cell <- if (identical(input$cell_mode, "file")) input$up_cell$datapath %||% "" else ""
      umap <- if (identical(input$umap_mode, "file")) input$up_umap$datapath %||% "" else ""
      open_dataset(u$objs, sub("^_|_$", "", stem), NULL, cell, umap)
    } else open_folder()
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
    for (e in list(figs, tables, gs_env, score_env)) rm(list = ls(e), envir = e)
    rv$figs_seen <- character(0); rv$tables_seen <- character(0)
    cache_env$data <- b$cache; b$cache <- NULL
    rm(list = ls(km_env), envir = km_env)
    D(b)
    if (nzchar(umap_file %||% "")) {
      xy <- tryCatch(read_umap(umap_file, b$mcs), error = function(e) { warn("UMAP file: ", conditionMessage(e), todo = "A UMAP is computed from the module scores instead when the UMAP page opens."); NULL })
      if (!is.null(xy)) UMAP(list(xy = xy, source = paste("From file:", basename(umap_file))))
    } else if (!is.null(cache_env$data$umap)) UMAP(cache_env$data$umap)
    if (nzchar(cell_file %||% "")) start_cells(cell_file)
  }

  start_cells <- function(path) {
    b <- D(); key <- paste0("cells:", file_key(path)); hit <- cache_env$data[[key]]
    if (!is.null(hit)) { CELLS(hit); rv$cell_status <- list(state = "done", file = basename(path)); return() }
    mcx_log("Reading cell information in the background: ", basename(path), " (", fmt_size(file.size(path)), ")")
    if (requireNamespace("callr", quietly = TRUE)) {
      f <- summarise_cells; environment(f) <- globalenv(); out <- tempfile(fileext = ".rds")
      job <- callr::r_bg(function(f, path, mcs, out) { saveRDS(f(path, mcs), out); TRUE }, args = list(f = f, path = path, mcs = b$mcs, out = out), supervise = TRUE)
      rv$job <- list(job = job, out = out, path = path, key = key)
      rv$cell_status <- list(state = "running", file = basename(path), size = file.size(path))
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
    rv$cell_status <- list(state = "done", file = basename(path))
    showNotification(HTML(paste0("<b>Cell information ready.</b> ", fmt_n(res$n_total), " cells. Labels and the Conditions page are now available.")), type = "message")
  }
  observe({
    j <- rv$job; req(j)
    if (j$job$is_alive()) { invalidateLater(1000); return() }
    res <- tryCatch({ j$job$get_result(); readRDS(j$out) }, error = function(e) e)
    rv$job <- NULL; finish_cells(res, j$key, j$path)
  })

  output$nav_status <- renderUI({
    b <- D(); s <- rv$cell_status
    if (is.null(b)) return(div(class = "mcx-navstatus", span(class = "mcx-chip muted", "No dataset open")))
    if (s$state == "running") return(div(class = "mcx-navstatus", span(class = "mcx-chip warn", icon("spinner", class = "fa-spin"), "Reading cell information")))
    div(class = "mcx-navstatus", span(class = paste("mcx-chip", if (s$state == "error") "warn"),
        title = switch(s$state, done = paste("Cell information:", s$file), error = paste("Cell information failed:", s$msg), "No cell information"),
        icon(switch(s$state, done = "circle-check", error = "triangle-exclamation", "layer-group")), b$name))
  })

  # ------------------------------------------------------------ HOME dashboard
  output$home_head <- renderUI({
    b <- req(D()); s <- rv$cell_status
    chip <- function(cls, ic, txt) span(class = paste("mcx-chip", cls), icon(ic), txt)
    div(class = "mcx-homehead",
      div(div(class = "mcx-eyebrow", "Home"), h2(class = "mcx-dsname", b$name), div(class = "mcx-path", icon("folder"), if (identical(b$src$counts, "(uploaded)")) "uploaded files" else dirname(b$src$counts))),
      div(class = "d-flex gap-2 flex-wrap",
        switch(s$state, done = chip("ok", "circle-check", "Cell information loaded"), running = chip("warn", "spinner", "Reading cell information"),
               error = chip("warn", "triangle-exclamation", "Cell information failed"), chip("muted", "circle-minus", "No cell information")),
        if (!is.null(UMAP())) chip("ok", "circle-check", "UMAP ready") else chip("muted", "braille", "UMAP made on first visit")))
  })
  output$home_tiles <- renderUI({
    b <- req(D()); cl <- CELLS(); pl <- primary_label()
    tile <- function(ic, v, l, sub = NULL) div(class = "mcx-tile", div(class = "mcx-tile-icon", icon(ic)),
      div(div(class = "v", v), div(class = "l", l), if (!is.null(sub)) div(class = "s", sub)))
    cond <- guess_cond(); unit <- if (!is.null(cond)) unit_candidates(cond)[1]; if (!length(unit) || is.na(unit)) unit <- NULL
    div(class = "mcx-tiles",
      tile("circle-nodes", fmt_n(length(b$mcs)), "Metacells", if (isTRUE(cl$cell_level)) paste0(formatC(cl$n_total / length(b$mcs), format = "f", digits = 1), " cells per metacell")),
      tile("circle", if (isTRUE(cl$cell_level)) fmt_n(cl$n_total) else "-", "Cells"),
      if (!is.null(unit)) tile("vial", ncol(cl$comp[[unit]]), pretty_name(unit)),
      if (!is.null(cond)) tile("vials", ncol(cl$comp[[cond]]), pretty_name(cond), paste(colnames(cl$comp[[cond]]), collapse = ", ")),
      if (!is.null(pl)) tile("tags", length(group_levels(pl)), key_label(pl)),
      tile("layer-group", ncol(b$Z), "Gene modules"),
      tile("table-cells", b$nb, "Pipeline blocks"),
      tile("dna", fmt_n(nrow(b$counts)), "Genes"))
  })
  output$home_cond <- renderPlot({
    cl <- req(CELLS()); cond <- req(guess_cond()); unit <- unit_candidates(cond)[1]; if (!length(unit) || is.na(unit)) unit <- NULL
    cells <- colSums(cl$comp[[cond]]); ns <- if (!is.null(unit)) { t <- cl$pairs[[unit]][[cond]]; t <- t[rowSums(t) > 0, , drop = FALSE]; table(factor(colnames(t)[max.col(t)], levels = names(cells))) } else NULL
    df <- data.frame(x = factor(names(cells), levels = names(cells)), y = as.vector(cells), n = if (!is.null(ns)) as.vector(ns[names(cells)]) else NA)
    ggplot(df, aes(x, y, fill = x)) + geom_col(width = 0.65) +
      geom_text(aes(label = paste0(fmt_n(y), if (!is.null(ns)) paste0("\n", n, " ", tolower(pretty_name(unit)), "s") else "")), vjust = -0.2, size = 3.3, lineheight = 0.9) +
      scale_fill_manual(values = pal_for(levels(df$x)), guide = "none") + scale_y_continuous(expand = expansion(mult = c(0, 0.28)), labels = fmt_n) +
      labs(x = NULL, y = "Cells") + theme_mcx(12)
  }, res = 96)
  output$home_lab_title <- renderUI({ pl <- primary_label(); if (is.null(pl)) "Cells per label" else paste("Cells per", key_label(pl)) })
  output$home_lab <- renderPlot({
    cl <- req(CELLS()); pl <- req(primary_label()); v <- mv(pl); req(!is.null(cl$n_cells))
    s <- tapply(cl$n_cells, v, sum); lv <- level_order(pl, names(s)); s <- s[lv]
    df <- data.frame(x = factor(names(s), levels = rev(names(s))), y = as.vector(s))
    ggplot(df, aes(y, x, fill = x)) + geom_col(width = 0.7) + scale_fill_manual(values = pal_for(levels(df$x)), guide = "none") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.08)), labels = fmt_n) + labs(x = "Cells", y = NULL) + theme_mcx(12) +
      theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "#ECEEF2"))
  }, res = 96)

  # ============================================================ when a dataset opens / cell info arrives
  observeEvent(D(), {
    b <- D()
    rv$ann <- data.frame(metacell = b$mcs, annotation = "undefined", stringsAsFactors = FALSE)
    rv$hist <- list(); rv$log <- data.frame(time = character(0), action = character(0), metacells = integer(0), name = character(0))
    rv$gops <- list(); rv$split_sig <- ""; rv$sel <- character(0); rv$labels_applied <- FALSE
    for (id in c("pl_mod", "cd_mod", "um_mod", "g_mod", "fm_mod")) set_mod_choices(id)
    set_mod_choices("h2_mods", b$mods$id[b$row_ord[1:min(4, ncol(b$Z))]]); set_mod_choices("hm_mods", b$mods$id[b$row_ord[1:min(6, ncol(b$Z))]])
    fg <- b$bs$gene_var; fg <- names(fg)[which.max(fg * (b$bs$gene_frac >= 5e-5))]
    for (id in c("pl_gene", "cd_gene", "fg_gene", "um_gene")) set_gene_choices(id, fg)
    updateNumericInput(session, "grp_k", value = min(max(b$nb, 2), 300), max = min(300, max(10, length(b$mcs) - 1)))
    refresh_choices()
    if (file.exists(autosave_file())) showNotification(
      tagList(paste0("Names saved automatically for this dataset on ", format(file.mtime(autosave_file()), "%d %b at %H:%M"), " were found. "),
              actionLink("ann_restore2", "Restore them")), duration = 20, type = "message", id = "autosave_note")
    rv$wizard <- FALSE
    opened$now <- TRUE; session$onFlushed(function() opened$now <- FALSE, once = TRUE)
  }, priority = 10)
  # when a dataset and its cached cell information arrive together, choices are refreshed once (a second pass would undo the defaults)
  opened <- new.env(); opened$now <- FALSE
  observeEvent(CELLS(), { if (!opened$now) refresh_choices() }, ignoreNULL = FALSE, ignoreInit = TRUE)

  refresh_choices <- function() {
    b <- D(); req(b); cl <- CELLS(); pl <- primary_label(); gcm <- group_choices(); mc <- meta_choices()
    first <- !is.null(pl) && !isTRUE(rv$labels_applied)
    if (first) rv$labels_applied <- TRUE
    keep_or <- function(id, def, valid = NULL) { cur <- isolate(input[[id]]); if (first || is.null(cur) || !nzchar(cur[1]) || (!is.null(valid) && !all(cur %in% valid))) def else cur }
    updateSelectInput(session, "grp_by", choices = c(list("Blocks" = c("Pipeline blocks" = "block", "New blocks (re-group the metacells)" = "cluster")),
                                                     if (length(mc)) list("Group by an existing label" = mc),
                                                     list("Other" = c("Names given so far" = "ann", "Every metacell on its own" = "none"))),
                      selected = isolate(input$grp_by) %||% "block")
    bars <- c(list("Annotate page" = c("Annotate blocks" = "grp", "Names given" = "ann"), "Pipeline" = c("Pipeline blocks" = "block")), if (length(mc)) list("Cell metadata" = mc))
    bar_def <- unique(c(pl, "ann"))
    updateSelectizeInput(session, "hm_bars", choices = bars, selected = if (first) bar_def else isolate(input$hm_bars) %||% "ann")
    updateSelectizeInput(session, "h2_bars", choices = bars, selected = if (first) unique(c("grp", pl)) else isolate(input$h2_bars) %||% "grp")
    updateSelectInput(session, "pl_x", choices = gcm, selected = keep_or("pl_x", pl %||% "grp"))
    updateSelectInput(session, "pl_split", choices = c(list("No split" = c("No split" = "none")), group_choices(sel = FALSE)), selected = isolate(input$pl_split) %||% "none")
    mch <- metric_choices(); updateSelectInput(session, "pl_metric", choices = mch, selected = keep_or("pl_metric", unname(if ("Cells" %in% names(mch)) mch["Cells"] else mch[1]), mch))
    updateSelectInput(session, "mg_x", choices = group_choices(), selected = keep_or("mg_x", pl %||% "grp"))
    updateSelectInput(session, "um_by", choices = c(gcm, list("Expression" = c("A module" = "module", "A gene" = "gene", "A signature" = "sig")), list("Metacell" = metric_choices())),
                      selected = keep_or("um_by", pl %||% "grp"))
    updateSelectInput(session, "tg_from", choices = group_choices(), selected = keep_or("tg_from", pl %||% "grp"))
    updateSelectInput(session, "vg_from", choices = group_choices(sel = FALSE), selected = keep_or("vg_from", "grp"))
    cc <- comp_cols()
    if (length(cc)) {
      updateSelectInput(session, "cd_cond", choices = setNames(cc, pretty_name(cc)), selected = keep_or("cd_cond", guess_cond(), cc))
      updateSelectInput(session, "cd_grp_by", choices = group_choices(), selected = keep_or("cd_grp_by", pl %||% "grp"))
    }
    disable("vg_go", is.null(cl$total_umis), "Needs UMI counts per cell (nCount_RNA) in the cell information, added on the Home page.")
  }

  # ============================================================ grouping engine
  # new blocks: the pipeline's own method (k-means as in ComplexHeatmap column_km), cached per number and seed
  km_env <- new.env()
  new_blocks <- function(k, seed) {
    key <- paste(k, seed); if (!is.null(km_env[[key]])) return(km_env[[key]])
    b <- D(); n <- nrow(b$Z)
    r <- withProgress(message = paste("Making", k, "new blocks"), detail = if (n > 5000) "a few seconds" else NULL, value = 0.5, km_blocks(b$Z, k, seed))
    assign(key, r, envir = km_env); r
  }

  grp_k_d <- debounce(reactive(input$grp_k), 500)
  observeEvent(input$grp_k_minus, updateNumericInput(session, "grp_k", value = max(2, (input$grp_k %||% 25) - 1)))
  observeEvent(input$grp_k_plus, updateNumericInput(session, "grp_k", value = min(300, (input$grp_k %||% 25) + 1)))
  base_groups <- reactive({
    b <- req(D()); by <- input$grp_by %||% "block"; ht_rank <- match(b$mcs, b$mc_order_ht)
    pending <- FALSE
    if (by == "cluster") {
      k <- suppressWarnings(as.integer(grp_k_d())); if (is.na(k)) k <- b$nb; k <- max(2, min(k, 300, length(b$mcs) - 1))
      seed <- suppressWarnings(as.integer(input$grp_seed)); if (is.na(seed)) seed <- 42L
      r <- new_blocks(k, seed); g <- r$id; ord <- r$order; first <- seq_len(max(g))
      out <- list(id = paste0("N", match(g, first)), levels = paste0("N", seq_along(first)), order = ord, sig = paste("cluster", k, seed))
    } else if (by == "none") {
      out <- list(id = b$mcs, levels = b$mc_order_ht, order = match(b$mc_order_ht, b$mcs), sig = "none")
    } else if (by == "block") {
      out <- list(id = paste0("B", b$block), levels = paste0("B", seq_len(b$nb)), order = match(b$mc_order_ht, b$mcs), sig = "block")
    } else {
      v <- mv(by); req(length(v) == length(b$mcs)); v <- as.character(v); v[is.na(v)] <- "NA"
      lv <- level_order(by, v)
      out <- list(id = v, levels = lv, order = order(match(v, lv), ht_rank), sig = paste("label", by))
    }
    out$by <- by; out$pending <- pending; out
  })
  grouping <- reactive({   # base blocks, then the splits and merges in the order they were made
    b <- D(); g <- base_groups(); ops <- if (identical(rv$split_sig, g$sig)) rv$gops else list()
    id <- g$id; lv <- g$levels; ord <- g$order; done <- character(0)
    for (op in ops) {
      if (op[["type"]] == "split") {
        G <- op[["g"]]; n <- op[["n"]]; if (!G %in% lv) next
        idx <- which(id == G); if (length(idx) <= n) next
        r <- km_blocks(b$Z[idx, , drop = FALSE], n, op[["seed"]] %||% 42L)   # same method as the blocks
        id[idx] <- paste0(G, ".", r$id)
        lv <- append(lv[lv != G], paste0(G, ".", seq_len(n)), after = match(G, lv) - 1)
        pos <- which(ord %in% idx); ord[pos] <- idx[r$order]
      } else {   # merge: members become one group at the place of the first; metacells keep their order inside
        gs <- lv[lv %in% op[["gs"]]]; if (length(gs) < 2) next
        id[id %in% gs] <- op[["name"]]
        lv <- append(lv[!lv %in% gs], op[["name"]], after = match(gs[1], lv) - 1)
        rk <- integer(length(ord)); rk[ord] <- seq_along(ord); ord <- order(match(id, lv), rk)
      }
      done <- c(done, paste(op[["type"]], op[["g"]] %||% "", op[["n"]] %||% "", paste(op[["gs"]], collapse = "+")))
    }
    names(id) <- b$mcs
    list(id = id, levels = lv, order = ord, sig = g$sig, by = g$by, mc_order = b$mcs[ord], nops = length(done),
         nsplit = sum(startsWith(done, "split")), nmerge = sum(startsWith(done, "merge")),
         key = paste(g$sig, paste(done, collapse = ","), length(lv), sum(match(id, lv) * (seq_along(id) %% 97))))
  })
  gstats_for <- function(id, lv, key) {
    b <- D()
    if (identical(key, "block")) return(b$bs)
    if (!is.null(gs_env[[key]])) return(gs_env[[key]])
    gs <- withProgress(message = "Finding the top genes of each group", value = 0.5, group_stats(b$counts, b$libsize, match(id, lv), length(lv)))
    if (length(ls(gs_env)) > 4) rm(list = ls(gs_env), envir = gs_env)
    assign(key, gs, envir = gs_env); gs
  }
  grouping_gstats <- function() { gt <- grouping(); gstats_for(gt$id, gt$levels, if (gt$sig == "block" && !gt$nops) "block" else gt$key) }
  grp_top_genes <- reactive({ gt <- grouping(); if (length(gt$levels) > 400) return(NULL); top_genes_all(grouping_gstats(), 10) })

  # ============================================================ ANNOTATE: heatmap
  observe({ disable_choice("hm_view", "grp", identical(input$grp_by, "none"),
                           "Blocks are switched off ('Every metacell on its own' is chosen in the Blocks bar), so there is nothing to average. Choose blocks there to use this.") })
  hm_base <- reactive({
    b <- req(D()); gt <- grouping(); view <- if (gt$by == "none") "mc" else input$hm_view %||% "mc"
    f <- factor(gt$id, gt$levels)
    if (view == "mc") { cols <- gt$mc_order; M <- t(b$Z[cols, , drop = FALSE]); grp <- unname(gt$id[cols]) }
    else { M <- t(rowsum(b$Z, f, reorder = TRUE) / tabulate(f, nlevels(f))); cols <- gt$levels; grp <- gt$levels }
    zs <- rv$zoom_set; zoom <- isTRUE(input$hm_zoom) && length(zs) > 0   # zoom: only the blocks that held the selection when zooming in
    if (zoom) { keep <- if (gt$by == "none") cols %in% zs else grp %in% unique(gt$id[zs])
      if (any(keep)) { cols <- cols[keep]; grp <- grp[keep]; M <- M[, keep, drop = FALSE] } else zoom <- FALSE }
    ro <- b$row_ord
    if (identical(input$hm_modsel, "pick")) { keep <- intersect(input$hm_mods, b$mods$id)
      validate(need(length(keep) > 0, "No modules chosen. Pick modules under 'Modules shown' on the left, or show all.")); ro <- ro[b$mods$id[ro] %in% keep] }
    M <- M[ro, , drop = FALSE]; rownames(M) <- b$mods$short[ro]
    bb <- label_bars(input$hm_bars, view, cols, f)
    list(M = M, cols = cols, grp = grp, bars = bb$bars, bc = bb$bc, view = view, gt = gt, ro = ro, zoom = zoom)
  })
  # genes of the modules shown, listed under the heatmap
  output$hm_modgenes <- renderDT({
    b <- req(D()); ro <- hm_base()$ro
    d <- data.frame(Module = b$mods$short[ro], `Heatmap row` = seq_along(ro), Genes = lengths(b$mod_genes[b$mods$id[ro]]),
                    `Genes in the module` = vapply(b$mod_genes[b$mods$id[ro]], paste, "", collapse = ", "), num = b$mods$num[ro], check.names = FALSE)
    datatable(d, rownames = FALSE, selection = "none", class = "compact stripe hover",
              options = list(dom = "ftp", pageLength = 10, order = list(list(1, "asc")),
                             columnDefs = list(list(targets = 3, className = "wrap"), list(targets = 0, orderData = 4), list(targets = 4, visible = FALSE, searchable = FALSE))))
  })
  output$hm_modgenes_n <- renderUI(paste0("Genes of the ", length(hm_base()$ro), " modules shown (search a gene or sort by module)"))
  # coloured label bars above a heatmap: one per chosen label
  label_bars <- function(keys, view, cols, f) {
    bars <- list(); bc <- list()
    for (k in keys) {
      v <- if (view == "mc") mv(k, cols) else group_major(k, f)[match(cols, levels(f))]
      if (is.null(v)) next
      v <- as.character(v); v[is.na(v)] <- "NA"; nm <- key_label(k)
      bars[[nm]] <- v; bc[[nm]] <- pal_for(level_order(k, v))
    }
    list(bars = bars, bc = bc)
  }
  legend_ui <- function(bc) {
    if (!length(bc)) return(NULL)
    tagList(lapply(names(bc), function(nm) {
      lv <- names(bc[[nm]]); show <- head(lv, 40)
      div(class = "mcx-legend", span(class = "t", nm),
          lapply(show, function(l) span(span(class = "swatch", style = paste0("background:", bc[[nm]][[l]])), l)),
          if (length(lv) > 40) span(paste0("and ", length(lv) - 40, " more")))
    }))
  }
  hm_opts <- reactive({ hc <- heat_cols()
    list(low = hc[1], mid = hc[2], high = hc[3], lim = iv("hm_lim", 2), sep = isTRUE(input$hm_sep), labels = isTRUE(input$hm_glab), height = iv("hm_height", 500)) })
  hm_draw <- function(W, H) {
    d <- hm_base(); o <- hm_opts(); none <- d$gt$by == "none"; o$width_px <- W; o$height <- H
    o$pr <- max(1, session$clientData$pixelratio %||% 1)
    heat_plot(d$M, if (none) rep("x", length(d$grp)) else d$grp, d$bars, d$bc, o)
  }
  hm_W <- reactive(session$clientData$output_hm_plot_width %||% 900)
  hm_plot <- reactive(hm_draw(hm_W(), iv("hm_height", 500)))
  # what the browser needs to turn the mouse into columns and blocks (hover box, click, drag), see assets/app.js
  hm_info <- reactive({
    d <- hm_base(); gt <- d$gt; b <- D(); none <- gt$by == "none"
    L <- heat_layout(c(rownames(d$M), names(d$bars)), if (none) rep("x", length(d$grp)) else d$grp, length(d$bars), hm_W(), iv("hm_height", 500), isTRUE(input$hm_glab))
    r <- rle(as.character(d$grp)); e <- cumsum(r$lengths); st <- e - r$lengths + 1
    f <- factor(gt$id, gt$levels); n <- tabulate(f, nlevels(f)); names(n) <- gt$levels
    pl <- primary_label()
    lab <- if (!none && !is.null(pl)) setNames(group_cell_major(pl, f), gt$levels)
    nm  <- if (!none) setNames(group_major("ann", f, pct = TRUE), gt$levels)
    info <- if (none) as.list(r$values) else lapply(r$values, function(g) paste0("<b>", g, "</b> · ", fmt_n(n[[g]]), " metacells",
      if (!is.null(lab)) paste0("<br>", htmltools::htmlEscape(key_label(pl)), ": ", htmltools::htmlEscape(lab[[g]])),
      paste0("<br>Name: ", if (is.na(nm[[g]]) || grepl(paste0("^(", paste(NO_LABEL, collapse = "|"), ")"), nm[[g]])) "<i>not named</i>" else htmltools::htmlEscape(nm[[g]]))))
    list(W = L$W, H = L$H, x0 = L$x0, x1 = L$x1, nc = L$nc, top = L$H - L$ytop, bottom = L$H - L$y0, zoom = isTRUE(d$zoom),
         runs = lapply(seq_along(st), function(i) list(s = st[i], e = e[i], info = info[[i]])))
  })
  observe({ session$sendCustomMessage("mcx-hm-geom", hm_info()) })
  # the selection is drawn by the browser on top of the image, so selecting never redraws the heatmap
  observe({
    d <- hm_base(); s <- rv$sel
    on <- if (!length(s)) logical(0) else if (d$view == "mc") d$cols %in% s else d$cols %in% unique(d$gt$id[s])
    rr <- rle(on); e <- cumsum(rr$lengths); st <- e - rr$lengths + 1
    session$sendCustomMessage("mcx-hm-sel", list(runs = lapply(which(rr$values), function(i) list(s = st[i], e = e[i]))))
  })
  output$hm_plot <- renderPlot({
    p <- hm_plot(); w <- input[["hmdl-w"]] %||% 12; h <- input[["hmdl-h"]] %||% 7
    reg_fig("01_heatmap", "Annotate", function() print(isolate(hm_draw(w * 96, h * 96))), w, h); p
  }, height = function() iv("hm_height", 500), res = 96)
  dl_server("hmdl", function() { w <- input[["hmdl-w"]] %||% 12; h <- input[["hmdl-h"]] %||% 7; p <- hm_draw(w * 96, h * 96); function() print(p) },
            function() paste0(D()$name, "_heatmap"))
  output$hm_legend <- renderUI(legend_ui(hm_base()$bc))
  observeEvent(input$hm_colors, { d <- hm_base()
    open_colors(unlist(lapply(d$bc, names)), "Heatmap colours", scale = list(ids = c("hm_low", "hm_mid", "hm_high"), def = c("#000080", "#FFFFFF", "#B22222"))) })
  hmo_draw <- reactive({
    b <- req(D()); ht <- b$ht
    ht <- tryCatch({ ht@ht_list[[1]]@column_title <- paste0("B", seq_len(b$nb)); ht }, error = function(e) b$ht)
    function() draw(ht)
  })
  output$hm_orig_plot <- renderPlot({ f <- hmo_draw(); f() }, height = function() iv("hm_height", 500), res = 96)
  dl_server("hmodl", hmo_draw, function() paste0(D()$name, "_pipeline_heatmap"))

  # ------------------------------------------------------------ blocks bar
  output$grp_note <- renderUI({
    req(D()); gt <- grouping(); by <- input$grp_by %||% "block"
    txt <- switch(by,
      block = "The blocks drawn by the pipeline (the column splits of its heatmap).",
      cluster = "Metacells re-grouped with the pipeline's own method (k-means on the module scores, seed below). The pipeline's number of blocks and seed 42 give the pipeline blocks back.",
      none = "No blocks: every metacell is its own column. Click or drag to select metacells.",
      ann = "Grouped by the names given so far.",
      paste0("Metacells grouped by ", key_label(by), " (the label carried by most of their cells)."))
    ch <- c(if (gt$nsplit) paste(gt$nsplit, if (gt$nsplit == 1) "split" else "splits"), if (gt$nmerge) paste(gt$nmerge, if (gt$nmerge == 1) "merge" else "merges"))
    tagList(paste0(fmt_n(length(gt$levels)), if (by == "none") " metacells. " else " blocks. "), txt,
            if (length(ch)) tagList(paste0(" ", paste(ch, collapse = " and "), ". "),
              with_tip(actionLink("gop_undo", "Undo last"), "Undo the last split or merge"), " · ",
              with_tip(actionLink("split_reset", "Undo all"), "Back to the blocks without any split or merge")))
  })

  # ------------------------------------------------------------ selection: click a block, drag a range; Shift/Ctrl adds
  group_mcs <- function(g) { gt <- grouping(); names(gt$id)[gt$id %in% g] }
  pick <- function(mcs, add) {
    if (isTRUE(add)) { rv$sel <- if (all(mcs %in% rv$sel)) setdiff(rv$sel, mcs) else union(rv$sel, mcs) }
    else rv$sel <- if (length(mcs) && setequal(mcs, rv$sel)) character(0) else mcs
  }
  observeEvent(input$hm_pick, {   # from the browser: columns c0..c1 under a click or a drag
    x <- input$hm_pick; d <- hm_base(); nc <- ncol(d$M)
    i0 <- max(1, min(x$c0, x$c1)); i1 <- min(nc, max(x$c0, x$c1)); if (i1 < i0) return()
    # with blocks, a click or drag takes every block it touches; without blocks, exactly the metacells under it
    mcs <- if (d$gt$by == "none") d$cols[i0:i1] else group_mcs(unique(d$grp[i0:i1]))
    if (i0 == i1) pick(mcs, x$add) else if (isTRUE(x$add)) rv$sel <- union(rv$sel, mcs) else rv$sel <- mcs
  })
  observeEvent(input$sel_esc, { rv$sel <- character(0) })
  observeEvent(input$hm_zoom, { rv$zoom_set <- if (isTRUE(input$hm_zoom)) rv$sel else character(0)
    if (isTRUE(input$hm_zoom) && !length(rv$sel)) note("Select one or more blocks first; the zoom then shows only them.") }, ignoreInit = TRUE)
  observeEvent(input$hm_dbl, {   # double-click: zoom into that block, or back out
    if (isTRUE(input$hm_zoom)) return(update_switch("hm_zoom", value = FALSE))
    d <- hm_base(); c <- input$hm_dbl$c; if (is.null(c) || c < 1 || c > ncol(d$M)) return()
    rv$sel <- if (d$gt$by == "none") d$cols[c] else group_mcs(d$grp[c]); update_switch("hm_zoom", value = TRUE)
  })
  observeEvent(input$sel_clear, { rv$sel <- character(0) })
  observeEvent(input$grp_all, { rv$sel <- D()$mcs })
  sel_groups <- reactive({   # blocks fully inside the selection
    gt <- grouping(); if (!length(rv$sel) || length(gt$levels) > 2000) return(character(0))
    tab <- table(factor(gt$id, gt$levels), names(gt$id) %in% rv$sel)
    if (!"TRUE" %in% colnames(tab)) return(character(0))
    rownames(tab)[tab[, "TRUE"] > 0 & (if ("FALSE" %in% colnames(tab)) tab[, "FALSE"] == 0 else TRUE)]
  })
  observeEvent(input$sel_more, {
    showModal(modalDialog(title = "Select metacells", size = "m", easyClose = TRUE,
      selectInput("sel_lab_col", "By a label", choices = group_choices(sel = FALSE), width = "100%"),
      selectizeInput("sel_lab_val", NULL, choices = NULL, multiple = TRUE, width = "100%", options = c(list(placeholder = "Values"), RM_BTN)),
      actionButton("sel_lab_go", "Select these", class = "btn-primary btn-sm"),
      hr(),
      textAreaInput("sel_names", "Or paste metacell names", rows = 3, width = "100%", placeholder = "Separated by commas, spaces or new lines"),
      actionButton("sel_names_go", "Select these", class = "btn-primary btn-sm"),
      footer = modalButton("Close")))
  })
  observeEvent(input$sel_lab_col, { req(D(), nzchar(input$sel_lab_col)); lv <- group_levels(input$sel_lab_col); updateSelectizeInput(session, "sel_lab_val", choices = lv, server = length(lv) > 300) })
  observeEvent(input$sel_lab_go, { req(length(input$sel_lab_val)); b <- D(); rv$sel <- b$mcs[mv(input$sel_lab_col) %in% input$sel_lab_val]; removeModal() })
  observeEvent(input$sel_names_go, {
    b <- D(); nm <- parse_list(input$sel_names); hit <- intersect(nm, b$mcs)
    rv$sel <- hit; removeModal()
    if (length(nm) > length(hit)) warn(length(nm) - length(hit), " names were not found in this dataset.", todo = "Metacell names are the column names of the counts file.")
  })

  # ------------------------------------------------------------ the selection panel (right of the heatmap)
  cell_shares <- function(key, mcs) {
    cl <- CELLS(); b <- D(); i <- match(mcs, b$mcs); type <- sub(":.*$", "", key); col <- sub("^[^:]*:", "", key)
    if (type == "maj" && !is.null(cl$comp[[col]])) x <- colSums(cl$comp[[col]][i, , drop = FALSE])
    else { v <- as.character(mv(key, mcs)); v[is.na(v)] <- "NA"; w <- if (!is.null(cl$n_cells)) cl$n_cells[i] else rep(1, length(i)); x <- tapply(w, v, sum) }
    x <- x[x > 0]; if (!length(x)) return(NULL); sort(x / sum(x), decreasing = TRUE)
  }
  sel_top_genes <- reactive({
    s <- rv$sel; req(length(s)); b <- D(); gt <- grouping(); g <- sel_groups()
    if (length(g) && gt$by != "none" && setequal(group_mcs(g), s)) {
      gs <- grouping_gstats(); r <- rank_genes(gs, match(g, gt$levels), min_frac = 5e-5)
    } else {
      if (length(s) == length(b$mcs)) return(character(0))
      gs <- withProgress(message = "Top genes of the selection", value = 0.5, group_stats(b$counts, b$libsize, ifelse(b$mcs %in% s, 1L, 2L), 2))
      r <- rank_genes(gs, 1, 2, min_frac = 5e-5)
    }
    head(r$gene, 8)
  })
  share_ui <- function(sh, n = 3) div(class = "mcx-share", lapply(head(seq_along(sh), n), function(k) tagList(
    span(class = "nm", names(sh)[k]), span(class = "pc", if (sh[k] < 0.005) "<1%" else paste0(round(100 * sh[k]), "%")),
    div(class = "bar", tags$i(style = sprintf("width:%.0f%%;background:%s", 100 * sh[k], pal_for(names(sh)[k])))))))
  insp_parts <- reactive({
    b <- req(D()); s <- rv$sel; gt <- grouping(); cl <- CELLS()
    if (!length(s)) return(list(head = div(class = "mcx-insp-empty", "Nothing selected yet.",
      tags$ol(tags$li("Click a block on the heatmap."), tags$li("Type a name below and press Enter."), tags$li("What the block contains is shown under the name box.")),
      "Drag across the heatmap to take several blocks at once.")))
    g <- sel_groups(); whole <- length(g) && gt$by != "none" && setequal(group_mcs(g), s)
    title <- if (whole) { if (length(g) == 1) g else paste0(length(g), " blocks: ", paste(head(g, 4), collapse = ", "), if (length(g) > 4) ", ...") }
             else paste(fmt_n(length(s)), if (length(s) == 1) "metacell" else "metacells")
    ncell <- if (!is.null(cl$n_cells)) sum(cl$n_cells[match(s, b$mcs)]) else NA
    sub <- paste0(if (whole) paste0(fmt_n(length(s)), " metacells") else "",
                  if (!is.na(ncell)) paste0(if (whole) " · " else "", fmt_n(ncell), " cells (", round(100 * ncell / sum(cl$n_cells), 1), "%)") else "")
    nm <- rv$ann$annotation[match(s, rv$ann$metacell)]; nmt <- sort(table(nm), decreasing = TRUE)
    named <- if (all(nm %in% NO_LABEL)) "not named yet" else if (length(nmt) == 1) names(nmt) else
      paste0(paste0(head(names(nmt), 3), " (", head(as.integer(nmt), 3), ")", collapse = ", "), if (length(nmt) > 3) ", ...")
    pl <- primary_label(); sh <- if (!is.null(pl)) cell_shares(pl, s)
    flag <- if (!is.null(sh)) { if (sh[1] >= 0.8) span(class = "mcx-flag clean", "one main label") else if (sh[1] < 0.6) span(class = "mcx-flag mixed", "mixed") }
    zm <- colMeans(b$Z[match(s, b$mcs), , drop = FALSE]); top <- order(-zm)[1:min(3, length(zm))]
    mods <- paste0(b$mods$short[top], " (", vapply(top, function(k) paste(head(b$mod_genes[[b$mods$id[k]]], 3), collapse = " "), ""), ")")
    tg <- tryCatch(sel_top_genes(), error = function(e) character(0))
    list(head = tagList(div(class = "mcx-insp-title", title, flag), div(class = "mcx-insp-sub", sub),
                        if (all(nm %in% NO_LABEL)) div(class = "mcx-kv mb-0", div(class = "k", "Names given"), div(class = "v", "not named yet"))
                        else { ns <- cell_shares("ann", s)
                          tagList(div(class = "mcx-kv mb-0", div(class = "k", "Names already given (share of cells)"), share_ui(ns, 4)),
                                  div(class = "mcx-why", icon("triangle-exclamation"), span(paste0(if (any(nm %in% NO_LABEL)) "Part of this selection is" else "This selection is",
                                    " already named. Applying a name replaces the names above for these metacells; Undo brings them back.")))) }),
         body = tagList(
      if (!is.null(sh)) div(class = "mcx-kv", div(class = "k", paste(key_label(pl), "(share of cells)")), share_ui(sh),
        { mt <- sort(table(as.character(mv(pl, s))), decreasing = TRUE)
          div(class = "mcx-mcnote", tags$b("Metacells by their own main label: "),
              paste0(paste0(fmt_n(head(as.integer(mt), 4)), " ", head(names(mt), 4), collapse = ", "), if (length(mt) > 4) " and others", ". The label bar colours each metacell by its own main label, so a block of mixed metacells can look uniform there while its cells are split.")) }),
      div(class = "mcx-kv", div(class = "k", "Highest modules"), div(class = "v", paste(mods, collapse = ", "))),
      if (length(tg)) div(class = "mcx-kv", div(class = "k", "Top genes (vs all other metacells)"), div(class = "v", paste(tg, collapse = ", "))),
      if (!is.null(sh) && sh[1] < 0.6 && whole && length(g) == 1) div(class = "mcx-why", icon("lightbulb"),
        span(paste0("Less than 60% of the cells share one ", key_label(pl), ". Splitting this block (below) may separate the cell types.")))))
  })
  output$insp <- renderUI(insp_parts()$head)
  output$insp_body <- renderUI(insp_parts()$body)
  output$can_split <- reactive({ gt <- grouping(); g <- sel_groups(); length(g) == 1 && gt$by != "none" && setequal(group_mcs(g), rv$sel) && sum(gt$id == g) > 2 })
  outputOptions(output, "can_split", suspendWhenHidden = FALSE)
  output$can_merge <- reactive({ gt <- grouping(); g <- sel_groups(); length(g) >= 2 && gt$by != "none" && setequal(group_mcs(g), rv$sel) })
  outputOptions(output, "can_merge", suspendWhenHidden = FALSE)
  output$merge_label <- renderUI({ g <- sel_groups(); paste0("Same population? Merge ", if (length(g) <= 4) paste(g, collapse = ", ") else paste(length(g), "blocks"), " into one") })
  output$split_label <- renderUI(paste0("Mixed? Split ", sel_groups()[1], " into"))

  # ------------------------------------------------------------ blocks table (rows mirror the selection)
  grp_proxy <- dataTableProxy("grp_tbl")
  grp_summary <- reactive({
    b <- req(D()); gt <- grouping(); cl <- CELLS(); lv <- gt$levels; f <- factor(gt$id, lv)
    n <- tabulate(f, length(lv)); pl <- primary_label()
    d <- data.frame(Block = lv, Metacells = n, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(cl$n_cells)) d$Cells <- as.vector(tapply(cl$n_cells, f, sum))
    if (gt$by == "none") return(d)
    if (any(!rv$ann$annotation %in% NO_LABEL) && gt$by != "ann") { x <- group_major("ann", f, pct = TRUE); d$`Name given` <- ifelse(grepl(paste0("^(", paste(NO_LABEL, collapse = "|"), ")( |$)"), x), "", x) }
    if (!is.null(pl) && gt$by != pl) d[[paste("Main", key_label(pl))]] <- group_major(pl, f, pct = TRUE)
    zm <- rowsum(b$Z, f, reorder = TRUE) / n
    d$`Highest modules` <- apply(zm, 1, function(r) paste(b$mods$short[match(colnames(zm)[order(-r)[1:min(3, length(r))]], b$mods$id)], collapse = ", "))
    tg <- grp_top_genes(); if (!is.null(tg)) d$`Top genes` <- vapply(tg, function(g) paste(head(g, 5), collapse = ", "), "")
    d
  })
  output$grp_tbl <- renderDT({
    d <- grp_summary(); rows <- which(d$Block %in% isolate(sel_groups())); rv$tbl_sync <- rows
    datatable(d, rownames = FALSE, selection = list(mode = "multiple", selected = rows), class = "compact stripe hover nowrap",
              options = list(dom = if (nrow(d) > 300) "ftp" else "t", pageLength = 300, paging = nrow(d) > 300, ordering = FALSE,
                             scrollY = "280px", scrollCollapse = TRUE, scrollX = TRUE, deferRender = TRUE))
  })
  observeEvent(sel_groups(), {
    d <- isolate(grp_summary()); rows <- which(d$Block %in% sel_groups())
    if (!identical(sort(rows), sort(rv$tbl_sync %||% integer(0)))) { rv$tbl_sync <- rows; selectRows(grp_proxy, rows) }
  }, ignoreNULL = FALSE)
  observeEvent(input$grp_tbl_rows_selected, {
    rows <- input$grp_tbl_rows_selected %||% integer(0); prev <- rv$tbl_sync %||% integer(0)
    if (identical(sort(rows), sort(prev))) return()
    lv <- grouping()$levels; added <- lv[setdiff(rows, prev)]; removed <- lv[setdiff(prev, rows)]
    rv$tbl_sync <- rows
    if (length(added)) rv$sel <- union(rv$sel, group_mcs(added))
    if (length(removed)) rv$sel <- setdiff(rv$sel, group_mcs(removed))
  }, ignoreNULL = FALSE)

  # ------------------------------------------------------------ splitting a mixed block
  observeEvent(input$split_go, {
    gt <- grouping(); g <- sel_groups()
    if (length(g) != 1) return(warn("Splitting works on one whole block.", todo = "Click a single block on the heatmap (without Shift or Ctrl), then Split."))
    n <- max(2, as.integer(input$split_n %||% 2))
    if (sum(gt$id == g) <= n) return(warn(g, " has only ", sum(gt$id == g), " metacells.", todo = "Choose fewer parts, or name it as it is."))
    if (!identical(rv$split_sig, gt$sig)) { rv$gops <- list(); rv$split_sig <- gt$sig }
    rv$gops <- c(rv$gops, list(list(type = "split", g = g, n = n, seed = suppressWarnings(as.integer(input$grp_seed)) %|na|% 42L)))
    rv$sel <- group_mcs(paste0(g, ".1"))
    note(g, " split into ", paste0(g, ".", seq_len(n), collapse = ", "), ". ", g, ".1 is selected.")
  })
  observeEvent(input$split_reset, { rv$gops <- list() })
  observeEvent(input$gop_undo, { rv$gops <- head(rv$gops, -1) })

  # ------------------------------------------------------------ merging blocks that belong together
  merge_name <- function(gs, lv) {
    nm <- if (length(gs) <= 3) paste(gs, collapse = "+") else paste0(gs[1], "+", length(gs) - 1, " more")
    while (nm %in% lv) nm <- paste0(nm, "'"); nm
  }
  observeEvent(input$merge_go, {
    gt <- grouping(); g <- sel_groups()
    if (gt$by == "none") return(warn("There are no blocks to merge.", todo = "Choose blocks in the Blocks bar first."))
    if (length(g) < 2 || !setequal(group_mcs(g), rv$sel)) return(warn("Merging works on two or more whole blocks.", todo = "Shift or Ctrl click the blocks to merge, then Merge."))
    g <- gt$levels[gt$levels %in% g]; nm <- merge_name(g, gt$levels)
    if (!identical(rv$split_sig, gt$sig)) { rv$gops <- list(); rv$split_sig <- gt$sig }
    rv$gops <- c(rv$gops, list(list(type = "merge", gs = g, name = nm)))
    rv$sel <- group_mcs(nm)
    note(paste(g, collapse = ", "), " merged into ", nm, ". Undo it with 'Undo last' in the Blocks bar.")
  })

  # ------------------------------------------------------------ naming
  autosave_file <- reactive(paste0(AUTOSAVE_DIR, "/", gsub("[^A-Za-z0-9_]+", "_", req(D())$name), "_annotation.csv"))
  log_it <- function(action, n, name) rv$log <- rbind(data.frame(time = format(Sys.time(), "%H:%M:%S"), action = action, metacells = as.integer(n), name = name), rv$log)
  set_ann <- function(new, action, n, name) {
    rv$hist <- c(list(rv$ann), head(rv$hist, 29)); rv$ann <- new; log_it(action, n, name)
    tryCatch({ dir.create(AUTOSAVE_DIR, recursive = TRUE, showWarnings = FALSE); utils::write.csv(new, autosave_file(), row.names = FALSE) },
             error = function(e) mcx_log("Autosave failed: ", conditionMessage(e)))
  }
  observe({ a <- rv$ann; req(a); session$sendCustomMessage("mcx-datalist", list(id = "ann_names", values = as.list(sort(setdiff(unique(a$annotation), NO_LABEL))))) })
  assign_name <- function(nm) {
    s <- rv$sel; nm <- trimws(nm %||% "")
    if (!length(s)) return(warn("Nothing is selected.", todo = "Click a block on the heatmap first; it is outlined and described in the Selection panel."))
    if (!nzchar(nm)) return(warn("The name box is empty.", todo = "Type a name, e.g. a cell type, then press Enter."))
    a <- rv$ann; a$annotation[a$metacell %in% s] <- nm; set_ann(a, "named", length(s), nm)
    rv$sel <- character(0); session$sendCustomMessage("mcx-clear", "ann_name")
    showNotification(HTML(paste0("<b>", fmt_n(length(s)), "</b> metacells named <b>", htmltools::htmlEscape(nm), "</b>. Click the next block.")), type = "message", duration = 3)
  }
  observeEvent(input$ann_assign, assign_name(input$ann_name))
  observeEvent(input$ann_enter, assign_name(input$ann_enter$name))
  observeEvent(input$ann_undo, {
    if (!length(rv$hist)) return(warn("Nothing to undo.", todo = "Undo reverts namings, renamings and removals made in this session."))
    rv$ann <- rv$hist[[1]]; rv$hist <- rv$hist[-1]; log_it("undo", NA, "")
  })
  ann_counts <- reactive({
    b <- req(D()); a <- req(rv$ann); cl <- CELLS()
    lv <- level_order("ann", a$annotation); f <- factor(a$annotation, lv)
    t <- data.frame(Name = lv, Metacells = tabulate(f, length(lv)), check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(cl$n_cells)) { t$Cells <- as.integer(tapply(cl$n_cells, f, sum)); t$`% of cells` <- round(100 * t$Cells / sum(t$Cells), 1) }
    t
  })
  output$ann_list <- renderUI({
    t <- ann_counts(); real <- t[!t$Name %in% NO_LABEL, , drop = FALSE]; und <- t[t$Name %in% NO_LABEL, , drop = FALSE]
    row <- function(nm, n, cells, undef = FALSE) div(class = paste("mcx-name-row", if (undef) "undef"), `data-name` = nm,
      span(class = "swatch", style = paste0("background:", pal_for(nm))),
      span(class = "nm", `data-act` = "pick", title = "Select these metacells", nm),
      span(class = "ct", paste0(fmt_n(n), " mc", if (!is.na(cells)) paste0(" · ", cells, "%"))),
      if (!undef) tagList(tags$button(type = "button", class = "act", `data-act` = "rename", title = "Rename", icon("pen")),
                          tags$button(type = "button", class = "act", `data-act` = "drop", title = "Remove this name (back to undefined)", icon("xmark"))))
    pc <- function(d) if (is.null(d$`% of cells`)) rep(NA, nrow(d)) else d$`% of cells`
    div(class = "mcx-names",
      if (!nrow(real)) div(class = "mcx-insp-empty", "No names yet. Named blocks are listed here."),
      lapply(seq_len(nrow(real)), function(i) row(real$Name[i], real$Metacells[i], pc(real)[i])),
      lapply(seq_len(nrow(und)), function(i) row(und$Name[i], und$Metacells[i], pc(und)[i], TRUE)))
  })
  observeEvent(input$ann_act, {
    x <- input$ann_act; nm <- x$name; req(nm)
    if (x$act == "pick") { rv$sel <- rv$ann$metacell[rv$ann$annotation == nm]; return() }
    if (x$act == "drop") { a <- rv$ann; i <- a$annotation == nm; a$annotation[i] <- "undefined"; set_ann(a, "removed", sum(i), nm)
      return(note("'", nm, "' removed; its ", sum(i), " metacells are undefined again. Undo brings it back.")) }
    if (x$act == "rename") { rv$ren <- nm
      showModal(modalDialog(title = paste("Rename", nm), easyClose = TRUE, textInput("ann_newname", "New name", value = nm, width = "100%"),
        footer = tagList(modalButton("Cancel"), actionButton("ann_rename_go", "Rename", class = "btn-primary")))) }
  })
  observeEvent(input$ann_rename_go, { nm <- rv$ren; new <- trimws(input$ann_newname %||% ""); req(nm, nzchar(new))
    a <- rv$ann; i <- a$annotation == nm; a$annotation[i] <- new
    if (!is.na(rv$colors[nm])) { rv$colors[new] <- rv$colors[nm] }
    set_ann(a, paste("renamed", nm), sum(i), new); removeModal() })
  observeEvent(input$ann_go_export, { updateRadioButtons(session, "ex_mode", selected = "ann"); go("export") })
  observeEvent(input$ann_start, {
    ch <- c(meta_choices(), "Pipeline blocks" = "block", "Annotate blocks" = "grp")
    showModal(modalDialog(title = "Start from an existing label", size = "m", easyClose = TRUE,
      p(class = "mcx-hint", "Copies a label (e.g. broad_label) as a first draft of the names, to be corrected block by block. This replaces the current names; Undo brings them back."),
      selectInput("ann_init_col", NULL, choices = ch, width = "100%"),
      footer = tagList(modalButton("Cancel"), actionButton("ann_init", "Copy", class = "btn-primary"))))
  })
  observeEvent(input$ann_init, {
    req(input$ann_init_col); a <- rv$ann
    v <- as.character(mv(input$ann_init_col, a$metacell)); v[is.na(v)] <- "undefined"; a$annotation <- v
    set_ann(a, paste("copied", key_label(input$ann_init_col)), nrow(a), ""); removeModal()
  })
  observeEvent(input$ann_load, {
    rv$ann_msg <- NULL
    showModal(modalDialog(title = "Load a saved annotation", size = "m", easyClose = TRUE,
      p(class = "mcx-hint", "The Excel file from this app, or a CSV with a column of metacell names and a column called 'annotation'."),
      fileInput("ann_file", NULL, accept = c(".xlsx", ".csv"), width = "100%"), uiOutput("ann_load_msg"),
      footer = modalButton("Close")))
  })
  load_ann <- function(path, label) {
    x <- tryCatch(read_annotation_file(path, D()$mcs), error = function(e) e)
    if (inherits(x, "error")) { rv$ann_msg <- conditionMessage(x); return() }
    a <- rv$ann; i <- match(a$metacell, x$metacell); a$annotation[!is.na(i)] <- x$annotation[i[!is.na(i)]]
    set_ann(a, paste("loaded", label), sum(!is.na(i)), "")
    rv$ann_msg <- paste0(fmt_n(sum(!is.na(i))), " metacells matched", if (any(!x$metacell %in% a$metacell)) paste0("; ", sum(!x$metacell %in% a$metacell), " names are not in this dataset") else "", ".")
  }
  observeEvent(input$ann_file, load_ann(input$ann_file$datapath, input$ann_file$name))
  output$ann_load_msg <- renderUI({ m <- rv$ann_msg; if (is.null(m)) NULL else div(class = "alert alert-info py-2 small", m) })
  restore_autosave <- function() { load_ann(autosave_file(), "autosave"); removeNotification("autosave_note"); note(rv$ann_msg %||% "Restored.") }
  observeEvent(input$ann_restore, restore_autosave())
  observeEvent(input$ann_restore2, restore_autosave())
  output$autosave_ui <- renderUI({
    f <- autosave_file()
    div(class = "mcx-hint mt-3", "Every change is also saved automatically to ", code(f), if (file.exists(f)) tagList(". ", actionLink("ann_restore", "Restore the last autosave")))
  })
  output$ann_log <- renderDT(datatable(rv$log %||% data.frame(), rownames = FALSE, selection = "none", class = "compact stripe",
                                       colnames = c("Time", "Change", "Metacells", "Name"), options = list(dom = "tp", pageLength = 10, ordering = FALSE)))
  observeEvent(input$ann_colors, open_colors(setdiff(level_order("ann", rv$ann$annotation), NO_LABEL), "Colours of the names"))

  # ------------------------------------------------------------ selection elsewhere
  sel_or_warn <- function() { if (!length(rv$sel)) { warn("Nothing is selected on the heatmap yet.", todo = "Click a block first, then use this link."); FALSE } else TRUE }
  observeEvent(input$go_sel_plot, { req(sel_or_warn()); g <- sel_groups()
    if (length(g) && grouping()$by != "none") set_groups("pl_x", "grp", "pl_groups", g) else set_groups("pl_x", "sel", "pl_groups", c("Selected", "Not selected"))
    go("plot") })
  observeEvent(input$go_sel_cond, { req(sel_or_warn()); if (is.null(CELLS())) return(warn("Conditions need cell information.", todo = "Add the Seurat object or a table of cells on the Home page."))
    g <- sel_groups(); updateRadioButtons(session, "cd_mode", selected = "prop")
    if (length(g) && grouping()$by != "none") set_groups("cd_grp_by", "grp", "cd_groups", g) else set_groups("cd_grp_by", "sel", "cd_groups", "Selected")
    go("cond") })
  observeEvent(input$go_sel_mixed, { g <- sel_groups(); req(length(g))
    updateRadioButtons(session, "mg_view", selected = "mixed")
    if (identical(input$vg_from, "grp")) updateSelectizeInput(session, "vg_groups", selected = g) else { rv$pending$vg_groups <- g; updateSelectInput(session, "vg_from", selected = "grp") }
    go("genes"); note("Click 'Find varying genes' to check ", paste(g, collapse = ", "), ".") })
  observeEvent(input$go_sel_markers, { req(sel_or_warn()); updateRadioButtons(session, "mg_view", selected = "markers")
    rv$pending$tg_a <- "Selected"; if (identical(input$tg_from, "sel")) updateSelectizeInput(session, "tg_a", selected = "Selected") else updateSelectInput(session, "tg_from", selected = "sel")
    go("genes") })

  # Excel: 'Annotation' for reading and plotting, 'Metacells' for mapping back onto cells
  ann_sheets <- function(top = TRUE) {
    b <- D(); cl <- CELLS(); a <- rv$ann; gt <- grouping()
    lv <- level_order("ann", a$annotation); f <- factor(a$annotation, lv); n <- tabulate(f, length(lv))
    s1 <- data.frame(Annotation = lv, Colour = unname(pal_for(lv)), Metacells = n, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(cl$n_cells)) { cells <- as.vector(tapply(cl$n_cells, f, sum)); s1$Cells <- cells; s1$`% of cells` <- round(100 * cells / sum(cells), 2) }
    pl <- primary_label()
    for (k in unique(c(pl, if (length(cl$labels)) paste0("lab:", names(cl$labels))))) s1[[paste("Main", key_label(k))]] <- group_cell_major(k, f)
    zm <- rowsum(b$Z, f, reorder = TRUE) / n
    s1$`Top modules` <- apply(zm, 1, function(r) { i <- match(colnames(zm)[order(-r)[1:min(3, length(r))]], b$mods$id)
      paste0(b$mods$short[i], " (", vapply(b$mod_genes[b$mods$id[i]], function(g) paste(head(g, 3), collapse = ", "), ""), ")", collapse = "; ") })
    if (top) { gs <- withProgress(message = "Preparing the Excel file", value = 0.5, group_stats(b$counts, b$libsize, as.integer(f), length(lv)))
               s1$`Top genes` <- vapply(top_genes_all(gs, 10), paste, "", collapse = ", ") }
    cc <- guess_cond()
    if (!is.null(cc)) { agg <- rowsum(cl$comp[[cc]], f, reorder = TRUE); pct <- round(100 * sweep(agg, 2, pmax(colSums(agg), 1), "/"), 2)
      for (j in colnames(pct)) s1[[paste0("% of ", j, " cells")]] <- pct[, j] }
    mc_col <- if (!is.null(cl$mc_col) && cl$mc_col != "..mc") cl$mc_col else "metacell"
    s2 <- data.frame(a$metacell, a$annotation, check.names = FALSE, stringsAsFactors = FALSE); names(s2) <- c(mc_col, "annotation")
    for (lc in names(cl$labels)) s2[[lc]] <- cl$labels[[lc]]
    if (isTRUE(startsWith(pl %||% "", "maj:"))) { lc <- sub("^maj:", "", pl); s2[[paste0("main_", lc)]] <- cl$comp_major[[lc]] }
    s2$pipeline_block <- b$block
    if (!gt$by %in% c("block", "none")) s2$heatmap_group <- unname(gt$id)
    if (!is.null(cl$n_cells)) s2$n_cells <- cl$n_cells
    list(Annotation = s1, Metacells = s2)
  }
  xlsx_handler <- function() downloadHandler(function() paste0(D()$name, "_annotation_", Sys.Date(), ".xlsx"), function(file) write_xlsx_sheets(file, ann_sheets()))
  output$ann_xlsx <- xlsx_handler(); output$ann_xlsx_exp <- xlsx_handler()
  ex_sheets <- reactive({ req(identical(input$nav, "export"), identical(input$ex_mode, "ann")); ann_sheets(FALSE) })
  output$ex_s1 <- renderDT(datatable(ex_sheets()$Annotation[, setdiff(names(ex_sheets()$Annotation), "Colour")], rownames = FALSE, selection = "none", class = "compact stripe nowrap",
                                     options = list(dom = "tp", pageLength = 10, scrollX = TRUE)), server = FALSE)
  output$ex_s2 <- renderDT(datatable(ex_sheets()$Metacells, rownames = FALSE, selection = "none", class = "compact stripe nowrap",
                                     options = list(dom = "ftp", pageLength = 10, scrollX = TRUE)), server = FALSE)

  # ============================================================ signatures
  sig_ids <- c("pl_sig", "cd_sig", "um_sig", "g_sig")
  open_sig_modal <- function(target = NULL, genes = character(0), name = "", method = "zmean", edit = NULL) {
    rv$sig_target <- target; rv$sig_edit <- edit
    showModal(modalDialog(title = if (is.null(edit)) "New gene signature" else paste("Edit", edit), size = "m", easyClose = TRUE,
      textInput("sig_name", "Name", value = name, width = "100%", placeholder = "e.g. IFNG response"),
      textAreaInput("sig_genes", "Genes", value = paste(genes, collapse = ", "), rows = 5, width = "100%", placeholder = "IFNG, STAT1, CXCL9, CXCL10, GBP1 ..."),
      uiOutput("sig_check"),
      radio_tips("sig_method", "Score", SIG_METHODS,
        list(zmean = "Each gene is z-scored across metacells, then averaged, so one highly expressed gene does not dominate.",
             ctrl = "Average expression of the genes minus the average of random genes with similar expression (Tirosh et al. 2016, Science; Seurat AddModuleScore)."),
        selected = method),
      footer = tagList(if (!is.null(edit)) actionButton("sig_delete", "Delete", icon = icon("trash"), class = "btn-light me-auto"),
                       modalButton("Cancel"), actionButton("sig_save", "Save signature", class = "btn-primary"))))
  }
  output$sig_check <- renderUI({
    b <- req(D()); m <- match_genes(parse_list(input$sig_genes), rownames(b$counts))
    if (!length(m$found) && !length(m$missing)) return(NULL)
    hint(paste0(length(m$found), " genes found"), if (length(m$missing)) paste0(". Not in this dataset: ", paste(head(m$missing, 12), collapse = ", "), if (length(m$missing) > 12) " ..."))
  })
  for (id in sig_ids) local({
    id <- id
    observeEvent(input[[paste0(id, "_new")]], open_sig_modal(id))
    observeEvent(input[[paste0(id, "_lib")]], open_lib(id))
    observeEvent(input[[paste0(id, "_edit")]], {
      nm <- input[[id]][1]; s <- rv$sigs[[nm %||% ""]]
      if (is.null(s)) return(warn("No signature is chosen.", todo = "Pick one in the list to edit it, or click + to make one."))
      open_sig_modal(id, s$genes, nm, s$method, edit = nm)
    })
    output[[paste0(id, "_info")]] <- renderUI({
      nm <- input[[id]]; if (!length(nm)) return(if (!length(rv$sigs)) hint("A signature is a gene list (e.g. an IFNG response) scored in every metacell. + types a list; the book opens published gene sets. Signatures are shared by every page."))
      s <- rv$sigs[[nm[1]]]; if (is.null(s)) return(NULL)
      if (length(nm) > 1) return(hint(paste(length(nm), "signatures chosen. The scoring method of each is shown when it is chosen alone.")))
      tagList(
        div(class = "mcx-scoring", tags$label(class = "control-label", lab("Scored as", SIG_METHOD_TIP)),
            selectInput(paste0(id, "_method"), NULL, SIG_METHODS, selected = s$method, width = "100%")),
        div(class = "mcx-genetext", span(class = "lbl", paste0(length(s$genes), " genes:")),
            paste(head(s$genes, 30), collapse = ", "), if (length(s$genes) > 30) " ...",
            if (length(s$missing)) tagList(br(), span(class = "text-muted", paste0(length(s$missing), " not measured here: ", paste(head(s$missing, 8), collapse = ", "), if (length(s$missing) > 8) " ...")))))
    })
    observeEvent(input[[paste0(id, "_method")]], {
      nm <- input[[id]][1]; m <- input[[paste0(id, "_method")]]; sg <- rv$sigs
      if (!is.null(sg[[nm %||% ""]]) && !identical(sg[[nm]]$method, m)) { sg[[nm]]$method <- m; rv$sigs <- sg }
    }, ignoreInit = TRUE)
  })
  observe({
    s <- names(rv$sigs); ss <- isolate(rv$sig_select)
    for (id in sig_ids) {
      cur <- isolate(input[[id]]); multi <- id %in% c("pl_sig", "cd_sig")
      sel <- intersect(cur, s)
      if (!is.null(ss) && identical(ss$id, id)) sel <- if (multi) unique(c(sel, ss$name)) else ss$name
      if (!length(sel) && length(s)) sel <- s[length(s)]
      updateSelectizeInput(session, id, choices = s, selected = sel)
    }
    rv$sig_select <- NULL
  })
  observeEvent(input$sig_save, {
    b <- D(); m <- match_genes(parse_list(input$sig_genes), rownames(b$counts)); nm <- trimws(input$sig_name %||% "")
    if (length(m$found) < 2) return(warn("A signature needs at least two genes found in this dataset.", todo = "Check the spelling; genes not found are listed under the box."))
    if (!nzchar(nm)) nm <- paste("Signature", length(rv$sigs) + 1)
    s <- rv$sigs; if (!is.null(rv$sig_edit) && rv$sig_edit != nm) s[[rv$sig_edit]] <- NULL
    s[[nm]] <- list(genes = m$found, missing = m$missing, method = input$sig_method %||% "zmean")
    rv$sig_select <- list(id = rv$sig_target, name = nm); rv$sigs <- s
    tg <- rv$sig_target
    if (!is.null(tg)) { w <- sub("_sig$", "_what", tg); if (w %in% c("pl_what", "cd_what")) updateRadioButtons(session, w, selected = "sig")
      if (tg == "um_sig") updateSelectInput(session, "um_by", selected = "sig"); if (tg == "g_sig") updateSelectInput(session, "g_src", selected = "sig") }
    else note("Signature '", nm, "' saved. It can be chosen under Signature on every page.")
    removeModal()
  })
  # gene-set library: MSigDB collections (msigdbr) or a .gmt file; ticked sets become signatures
  gs_cache <- new.env()
  open_lib <- function(target) {
    rv$sig_target <- target; ok <- genesets_ok()
    showModal(modalDialog(title = "Gene-set library", size = "l", easyClose = TRUE,
      p(class = "mcx-hint", "Published gene sets, each added as a signature and scored in every metacell. Genes not measured in this dataset are left out of the score."),
      if (!ok) div(class = "mcx-why", icon("circle-info"), span("MSigDB collections need the msigdbr package. In R: install.packages(\"msigdbr\"), then restart the app. A .gmt file (from MSigDB, Enrichr or GSEA) works without it.")),
      layout_columns(col_widths = c(7, 5),
        selectInput("lib_src", "Collection", c(GENESET_COLLECTIONS, "A .gmt file" = "gmt"), selected = if (ok) "H" else "gmt", width = "100%"),
        selectInput("lib_species", lab("Species", "Guessed from the gene names (upper case: human)."), c("Human" = "Homo sapiens", "Mouse" = "Mus musculus"),
                    selected = guess_species(rownames(D()$counts)), width = "100%")),
      conditionalPanel("input.lib_src == 'gmt'", fileInput("lib_gmt", NULL, accept = c(".gmt", ".txt"), buttonLabel = "Choose a .gmt file", width = "100%")),
      uiOutput("lib_msg"), DTOutput("lib_tbl"),
      footer = tagList(div(class = "me-auto mcx-scoring-inline", span(class = "control-label", lab("Score the added sets as", SIG_METHOD_TIP)),
                           selectInput("lib_method", NULL, SIG_METHODS, width = "260px", selectize = FALSE)),
                       modalButton("Close"), actionButton("lib_add", "Add ticked sets as signatures", icon = icon("plus"), class = "btn-primary"))))
  }
  lib_sets <- reactive({
    src <- req(input$lib_src)
    if (src == "gmt") { f <- input$lib_gmt; if (is.null(f)) return(NULL); return(tryCatch(read_gmt(f$datapath), error = function(e) e)) }
    key <- paste(src, input$lib_species)
    if (is.null(gs_cache[[key]])) {
      x <- withProgress(message = "Loading gene sets", value = 0.5, tryCatch(load_genesets(src, input$lib_species), error = function(e) e))
      if (inherits(x, "error")) return(x); assign(key, x, envir = gs_cache)
    }
    gs_cache[[key]]
  })
  lib_df <- reactive({
    x <- lib_sets(); if (is.null(x) || inherits(x, "error")) return(NULL); u <- rownames(D()$counts)
    miss <- lapply(x, function(g) g[!g %in% u])
    data.frame(`Gene set` = nice_set_name(names(x)), Genes = lengths(x), `Measured here` = lengths(x) - lengths(miss),
               `Not measured here` = vapply(miss, function(g) if (!length(g)) "" else paste0(paste(head(g, 10), collapse = ", "), if (length(g) > 10) paste0(" +", length(g) - 10, " more")), ""),
               genes = vapply(x, paste, "", collapse = " "), check.names = FALSE)
  })
  output$lib_msg <- renderUI({
    x <- lib_sets()
    if (is.null(x)) return(hint("Choose a .gmt file: one gene set per line, name, description, then the genes, separated by tabs."))
    if (inherits(x, "error")) return(div(class = "mcx-why", icon("triangle-exclamation"), span(conditionMessage(x))))
    hint(fmt_n(length(x)), " gene sets. Search by a word of the name (interferon) or by genes (IFNG CXCL9): sets containing them are listed. Click rows to tick them, then add.")
  })
  output$lib_tbl <- renderDT({ d <- lib_df(); req(d)
    datatable(d, rownames = FALSE, selection = "multiple", class = "compact stripe hover",
              options = list(dom = "ftp", pageLength = 8, order = list(), search = list(smart = TRUE),
                             columnDefs = list(list(targets = 4, visible = FALSE, searchable = TRUE), list(targets = 3, className = "wrap small")))) })
  observeEvent(input$lib_add, {
    rows <- input$lib_tbl_rows_selected; x <- lib_sets()
    if (!length(rows) || is.null(x) || inherits(x, "error")) return(warn("No gene sets are ticked.", todo = "Click rows of the table to tick them, then add."))
    b <- D(); s <- rv$sigs; added <- character(0)
    for (r in head(rows, 20)) {
      m <- match_genes(x[[r]], rownames(b$counts)); if (length(m$found) < 2) next
      nm <- nice_set_name(names(x)[r]); s[[nm]] <- list(genes = m$found, missing = m$missing, method = input$lib_method %||% "zmean"); added <- c(added, nm)
    }
    if (!length(added)) return(warn("None of the ticked sets has two or more genes measured here.", todo = "Check the species, or choose other sets."))
    rv$sig_select <- list(id = rv$sig_target, name = added[length(added)]); rv$sigs <- s
    tg <- rv$sig_target
    if (!is.null(tg)) { w <- sub("_sig$", "_what", tg); if (w %in% c("pl_what", "cd_what")) updateRadioButtons(session, w, selected = "sig")
      if (tg == "um_sig") updateSelectInput(session, "um_by", selected = "sig"); if (tg == "g_sig") updateSelectInput(session, "g_src", selected = "sig") }
    removeModal(); note(length(added), if (length(added) == 1) " signature added: " else " signatures added: ", paste(head(added, 3), collapse = "; "), if (length(added) > 3) " ...")
  })
  observeEvent(input$sig_delete, { s <- rv$sigs; s[[rv$sig_edit]] <- NULL; rv$sigs <- s; removeModal() })
  sig_scores <- function(nm) {
    s <- rv$sigs[[nm]]; if (is.null(s)) return(NULL)
    key <- paste(s$method, paste(sort(s$genes), collapse = ","))
    if (is.null(score_env[[key]])) assign(key, score_signature(D(), s$genes, s$method), envir = score_env)
    score_env[[key]]
  }

  # ============================================================ measuring per metacell (one or more features)
  measure <- function(p, mcs) {
    b <- D(); what <- input[[paste0(p, "_what")]] %||% "module"
    if (what == "module") {
      ids <- input[[paste0(p, "_mod")]]; ids <- ids[ids %in% b$mods$id]
      validate(need(length(ids) > 0, "Choose a module: type its number or a gene it contains."))
      i <- match(ids, b$mods$id); V <- b$Z[mcs, ids, drop = FALSE]; colnames(V) <- b$mods$lab[i]
      return(list(V = V, pos = V > 0, title = if (length(ids) == 1) b$mods$title[i] else NULL, short = paste(b$mods$short[i], collapse = "_"),
                  ylab = "Module score (z)", pos_name = "% metacells\nabove average", kind = "module", items = ids))
    }
    if (what == "gene") {
      g <- input[[paste0(p, "_gene")]]; g <- g[g %in% rownames(b$counts)]; validate(need(length(g) > 0, "Type a gene name."))
      V <- t(as.matrix(log_expr(b, g, mcs))); colnames(V) <- g
      return(list(V = V, pos = V > 0, title = if (length(g) == 1) g else NULL, short = paste(head(g, 4), collapse = "_"),
                  ylab = "Expression, log2(1 + UMIs per 100k)", pos_name = "% metacells\ndetected", kind = "gene", items = g))
    }
    if (what == "sig") {
      nms <- input[[paste0(p, "_sig")]]; nms <- nms[nms %in% names(rv$sigs)]
      validate(need(length(nms) > 0, "Make a signature first: click New and paste the genes."))
      V <- vapply(nms, function(n) unname(sig_scores(n)[mcs]), numeric(length(mcs))); V <- matrix(V, ncol = length(nms), dimnames = list(mcs, nms))
      return(list(V = V, pos = V > 0, title = if (length(nms) == 1) paste0(nms, " (", length(rv$sigs[[nms]]$genes), " genes)") else NULL, short = safe_name(paste(nms, collapse = "_")),
                  ylab = "Signature score", pos_name = "% metacells\nabove average", kind = "sig", items = nms))
    }
    k <- input[[paste0(p, "_metric")]]; validate(need(nzchar(k %||% ""), "Choose a metric."))
    v <- mv(k, mcs); validate(need(!is.null(v), "This metric needs cell information."))
    V <- matrix(as.numeric(v), ncol = 1, dimnames = list(mcs, key_label(k)))
    list(V = V, pos = V > stats::median(V, na.rm = TRUE), title = key_label(k), short = safe_name(key_label(k)), ylab = key_label(k), pos_name = "% above median", kind = "metric", items = k)
  }
  long_df <- function(m, ids) data.frame(id = rep(ids, ncol(m$V)), feature = factor(rep(colnames(m$V), each = length(ids)), levels = colnames(m$V)),
                                        y = as.vector(m$V), pos = as.vector(m$pos), stringsAsFactors = FALSE)
  for (p in c("pl", "cd")) local({
    p <- p
    output[[paste0(p, "_modinfo")]] <- renderUI({
      b <- req(D()); ids <- input[[paste0(p, "_mod")]]; ids <- ids[ids %in% b$mods$id]; if (!length(ids)) return(NULL)
      if (length(ids) == 1) return(genes_ui(b$mod_genes[[ids]], 40))
      hint(paste(b$mods$short[match(ids, b$mods$id)], collapse = ", "), " are drawn side by side.")
    })
  })
  stat_opts <- function(p, split = FALSE) list(cmp = input[[paste0(p, "_cmp")]] %||% "all", ref = input[[paste0(p, "_ref")]], pairs = rv$pairs[[p]],
    test = input[[paste0(p, "_test")]] %||% "wilcox", adj = input[[paste0(p, "_adj")]] %||% "holm", label = input[[paste0(p, "_lab")]] %||% "p",
    hidens = isTRUE(input[[paste0(p, "_hidens")]]), size = iv(paste0(p, "_statsize"), 3.5), split = split)
  style_opts <- function(p) list(points = isTRUE(input[[paste0(p, "_points")]]), pt = iv(paste0(p, "_pt"), 1.2), alpha = iv(paste0(p, "_alpha"), 0.6),
    mean = isTRUE(input[[paste0(p, "_mean")]]), angle = input[[paste0(p, "_angle")]] %||% "auto", txt = iv(paste0(p, "_txt"), 12),
    legend = isTRUE(input[[paste0(p, "_legend")]]), values = isTRUE(input[[paste0(p, "_values")]]), scale = isTRUE(input[[paste0(p, "_scale")]]))
  stats_dt <- function(s) datatable(s, rownames = FALSE, selection = "none", class = "compact stripe", options = list(dom = "tp", pageLength = 8, scrollX = TRUE)) |>
    formatSignif(intersect(c("p", "Adjusted p"), names(s)), 3)
  row_tests <- function(df, test) {   # one test per row of a heatmap / dot plot
    sp <- split(df, df$feature, drop = TRUE)
    pv <- vapply(sp, function(d) p_overall(d$y, d$x, test), 0)
    data.frame(Row = names(sp), Test = if (nlevels(droplevels(df$x)) == 2) (if (test == "t") "t-test" else "Wilcoxon") else (if (test == "t") "ANOVA" else "Kruskal-Wallis"),
               Columns = nlevels(droplevels(df$x)), p = pv, `Adjusted p` = stats::p.adjust(pv, "BH"), Significance = sig_stars(stats::p.adjust(pv, "BH")), check.names = FALSE)
  }
  # comparison choices follow the groups being compared; chosen pairs are added with two menus
  stats_sync <- function(p, levels_fun) {
    observe({
      lv <- levels_fun(); cur_ref <- isolate(input[[paste0(p, "_ref")]])
      updateSelectInput(session, paste0(p, "_ref"), choices = lv, selected = if (isTRUE(cur_ref %in% lv)) cur_ref else lv[1])
      updateSelectInput(session, paste0(p, "_pa"), choices = lv, selected = lv[1]); updateSelectInput(session, paste0(p, "_pb"), choices = lv, selected = lv[min(2, length(lv))])
      if (length(lv) > 3 && isolate(input[[paste0(p, "_cmp")]]) == "all") { updateSelectInput(session, paste0(p, "_cmp"), selected = "ref"); update_switch(paste0(p, "_hidens"), value = TRUE) }
      pr <- rv$pairs[[p]]; keep <- pr[vapply(strsplit(pr, " vs ", fixed = TRUE), function(z) all(z %in% lv), TRUE)]
      if (!identical(keep, pr)) rv$pairs[[p]] <- keep
    })
    observeEvent(input[[paste0(p, "_padd")]], {
      a <- input[[paste0(p, "_pa")]]; z <- input[[paste0(p, "_pb")]]
      if (is.null(a) || is.null(z) || a == z) return(warn("A pair needs two different groups.", todo = "Pick another group in the second box."))
      pr <- paste(a, "vs", z); if (!pr %in% rv$pairs[[p]] && !paste(z, "vs", a) %in% rv$pairs[[p]]) rv$pairs[[p]] <- c(rv$pairs[[p]], pr)
    })
    observeEvent(input[[paste0(p, "_pdel")]], { rv$pairs[[p]] <- setdiff(rv$pairs[[p]], input[[paste0(p, "_pdel")]]) })
    output[[paste0(p, "_pairs_ui")]] <- renderUI({
      pr <- rv$pairs[[p]]; if (!length(pr)) return(hint("No pairs yet. Choose two groups above and click +."))
      div(class = "mcx-pairlist", lapply(pr, function(x) span(class = "mcx-pair", x,
        tags$a(href = "#", class = "mcx-pair-x", title = "Remove", onclick = sprintf("Shiny.setInputValue('%s_pdel', %s, {priority:'event'});return false;", p, jsonlite::toJSON(x, auto_unbox = TRUE)), icon("xmark")))))
    })
  }
  heights <- function(h, type, k, ncol = 4) {
    if (type %in% c("heatmap", "dotplot")) return(max(if (type == "dotplot") 330 else 280, min(h * 1.8, 190 + 30 * k)))
    if (k > ncol) h * ceiling(k / ncol) * 0.62 else h
  }

  # ============================================================ EXPRESSION BY GROUP
  pl_lv <- reactive({ req(D(), input$pl_x); list(key = input$pl_x, lv = group_levels(input$pl_x)) })
  observeEvent(pl_lv(), fill_groups_rank("pl_groups", pl_lv()$key))
  observeEvent(rv$colors, { if (!is.null(D()) && !is.null(input$pl_x)) fill_groups_rank("pl_groups", input$pl_x, ticked = input$pl_groups) }, ignoreInit = TRUE)
  pl_cmp_levels <- reactive({
    if ((input$pl_split %||% "none") == "none") return(input$pl_groups %||% character(0))
    group_levels(input$pl_split)
  })
  stats_sync("pl", pl_cmp_levels)
  # statistics settings appear (and open) only while Statistics is on
  for (p in c("pl", "cd")) local({ p <- p; first <- TRUE
    observeEvent(input[[paste0(p, "_stats")]], { on <- isTRUE(input[[paste0(p, "_stats")]])
      session$sendCustomMessage("mcx-acc", list(acc = paste0(p, "_acc"), panel = "stats", show = on, open = on && !first)); first <<- FALSE }, ignoreNULL = FALSE)
  })
  # metrics (cells, UMIs...) are one value per metacell: heatmap, dot plot and row scaling need modules, genes or signatures
  observe({
    met <- identical(input$pl_what, "metric"); why <- "Needs modules, genes or signatures; metrics such as Cells are shown as box, violin, bar or strip plots."
    for (v in c("heatmap", "dotplot")) disable_choice("pl_type", v, met, why)
    if (met && isTRUE(input$pl_type %in% c("heatmap", "dotplot"))) updateRadioButtons(session, "pl_type", selected = "box")
  })
  pl_nfeat <- reactive(switch(input$pl_what %||% "module", module = length(input$pl_mod), gene = length(input$pl_gene), sig = length(input$pl_sig), 1))

  pl_result <- reactive({
    b <- req(D()); xk <- req(input$pl_x); groups <- input$pl_groups
    validate(need(length(groups) > 0, "Tick at least one group in the list on the left."))
    x <- as.character(mv(xk)); validate(need(length(x) > 0, "These groups need cell information.")); x[is.na(x)] <- "NA"
    groups <- intersect(groups, unique(x)); validate(need(length(groups) > 0, "None of the ticked groups are in the data."))
    sel <- b$mcs[x %in% groups]
    m <- measure("pl", sel); df <- long_df(m, sel); df$x <- factor(x[match(df$id, b$mcs)], levels = groups)
    type <- input$pl_type %||% "box"; k <- ncol(m$V); disp <- disp_labels(xk, groups); so <- style_opts("pl"); hc <- heat_cols()
    nm <- paste0(b$name, "_", m$short, "_by_", safe_name(key_label(xk)))
    if (type %in% c("heatmap", "dotplot")) {
      ref <- if (isTRUE(input$pl_ref %in% groups)) input$pl_ref else groups[1]
      st <- if (isTRUE(input$pl_stats) && length(groups) > 1) tile_tests(df, ref, input$pl_test %||% "wilcox", input$pl_adj %||% "holm", input$pl_lab, isTRUE(input$pl_hidens)) else NULL
      p <- summary_plot(df, type, c(so, list(size_by = "pos", size_name = m$pos_name, val_name = m$ylab, xlabels = disp, low = hc[1], mid = hc[2], high = hc[3],
                                            title = m$title, tile_labels = st, ref = if (!is.null(st)) ref, statsize = iv("pl_statsize", 3.5))))
      return(list(plot = p, data = df, stats = if (!is.null(st)) st[, setdiff(names(st), "label")], name = nm, k = k))
    }
    split <- (input$pl_split %||% "none") != "none" && input$pl_split != xk
    if (split) { s <- as.character(mv(input$pl_split, df$id)); df$fill <- factor(s, levels = level_order(input$pl_split, s)); df <- df[!is.na(df$fill), ] }
    else df$fill <- df$x
    if (k > 1) df$facet <- df$feature
    so$split <- split; so$xlabels <- disp; if (!split) so$fill_labels <- disp
    fl <- key_label(if (split) input$pl_split else xk)
    p <- dist_plot(df, type, pal_for(levels(df$fill)), so) + labs(title = m$title, x = NULL, y = m$ylab, fill = fl, colour = fl)
    if (k > 1) p <- p + facet_wrap(~ facet, scales = "free_y", ncol = min(4, k))
    st <- NULL
    if (isTRUE(input$pl_stats)) {
      o <- stat_opts("pl", split); st <- run_stats(df, o); p <- add_brackets(p, st, o, if (k > 1) levels(df$feature))
      if (!split && k == 1 && !is.null(st$overall)) p <- p + labs(subtitle = paste0(st$overall$Test[1], " across groups: p = ", fmt_p(st$overall$p[1])))
      st <- stats_table(st)
    }
    list(plot = p, data = df, stats = st, name = nm, k = k)
  })
  output$pl_plot <- renderPlot({
    r <- pl_result(); reg_fig("03_expression_by_group", "Expression by group", function() print(r$plot), input[["pldl-w"]], input[["pldl-h"]])
    if (!is.null(r$stats)) reg_table("expression_by_group_statistics", r$stats)
    r$plot
  }, height = function() heights(iv("pl_height", 480), input$pl_type %||% "box", pl_nfeat()), res = 96)
  dl_server("pldl", function() { p <- pl_result()$plot; function() print(p) }, function() pl_result()$name)
  output$pl_stats_tbl <- renderDT({ s <- pl_result()$stats; validate(need(!is.null(s), "Needs at least two groups with two or more metacells.")); stats_dt(s) })
  output$pl_dl_stats <- downloadHandler(function() paste0(pl_result()$name, "_stats.xlsx"), function(file) write_xlsx_sheets(file, list(Statistics = pl_result()$stats %||% data.frame())))
  observeEvent(input$pl_colors, { d <- pl_result()$data; open_colors(if (!is.null(d$fill)) levels(d$fill) else levels(d$x)) })

  # same features elsewhere
  copy_features <- function(from, to) {
    what <- input[[paste0(from, "_what")]] %||% "module"; if (what == "metric") return(FALSE)
    updateRadioButtons(session, paste0(to, "_what"), selected = what)
    if (what == "module") set_mod_choices(paste0(to, "_mod"), input[[paste0(from, "_mod")]])
    if (what == "gene") set_gene_choices(paste0(to, "_gene"), input[[paste0(from, "_gene")]])
    if (what == "sig") updateSelectizeInput(session, paste0(to, "_sig"), selected = input[[paste0(from, "_sig")]])
    TRUE
  }
  observeEvent(input$pl_go_cond, {
    if (is.null(CELLS())) return(warn("Conditions need cell information.", todo = "Add the Seurat object or a table of cells on the Home page."))
    if (!copy_features("pl", "cd")) return(warn("Metrics are not available on the Conditions page.", todo = "Choose a module, gene or signature under What to plot."))
    updateRadioButtons(session, "cd_mode", selected = "expr")
    xk <- input$pl_x; g <- input$pl_groups
    if (!identical(xk, "sel") && length(g) && length(g) < length(group_levels(xk))) { rv$pending$cd_cells <- g; if (identical(input$cd_grp_by, xk)) updateSelectizeInput(session, "cd_cells", selected = g) else updateSelectInput(session, "cd_grp_by", selected = xk) }
    go("cond")
  })
  observeEvent(input$pl_go_umap, {
    what <- input$pl_what
    if (what == "module") { updateSelectInput(session, "um_by", selected = "module"); set_mod_choices("um_mod", input$pl_mod[1]) }
    else if (what == "gene") { updateSelectInput(session, "um_by", selected = "gene"); set_gene_choices("um_gene", input$pl_gene[1]) }
    else if (what == "sig") { updateSelectInput(session, "um_by", selected = "sig"); updateSelectizeInput(session, "um_sig", selected = input$pl_sig[1]) }
    else updateSelectInput(session, "um_by", selected = input$pl_metric)
    go("umap")
  })

  # ============================================================ CONDITIONS
  observeEvent(input$cd_cond, {
    cl <- req(CELLS()); cond <- input$cd_cond; req(cond %in% comp_cols())
    lv <- colnames(cl$comp[[cond]]); rank_set("cd_conds", lv, lv, TRUE, pal_for(lv))
    uc <- unit_candidates(cond); updateSelectInput(session, "cd_unit", choices = setNames(uc, pretty_name(uc)), selected = uc[1])
    fc <- setdiff(comp_cols(), cond); updateSelectInput(session, "cd_filt_col", choices = c("(all samples)" = "none", setNames(fc, pretty_name(fc))))
  })
  observeEvent(input$cd_filt_col, {
    cl <- req(CELLS()); fc <- input$cd_filt_col
    updateSelectizeInput(session, "cd_filt_val", choices = if (fc %in% comp_cols()) colnames(cl$comp[[fc]]) else character(0))
  })
  cd_grp_levels <- reactive({ req(D(), input$cd_grp_by); group_levels(input$cd_grp_by) })
  observeEvent(cd_grp_levels(), {
    lv <- cd_grp_levels(); dl <- disp_labels(input$cd_grp_by, lv); real <- setdiff(lv, NO_LABEL)
    fill_groups_rank("cd_groups", input$cd_grp_by, default_n = 2)
    pend <- isolate(rv$pending$cd_cells); rv$pending$cd_cells <- NULL
    updateSelectizeInput(session, "cd_cells", choices = setNames(lv, dl[lv]), selected = pend %||% real[1], server = length(lv) > 300)
  })
  stats_sync("cd", reactive(if (identical(input$cd_type, "stacked")) character(0) else input$cd_conds %||% character(0)))
  observe({ ex <- identical(input$cd_mode, "expr"); disable_choice("cd_type", "stacked", ex, "Stacked bars show cell composition. Choose Cell composition at the top to use them.")
    if (ex && identical(isolate(input$cd_type), "stacked")) updateRadioButtons(session, "cd_type", selected = "box") })
  observe({ disable("cd_stats", identical(input$cd_type, "stacked"), "Stacked bars have no statistics. Choose Box, Violin, Bar or Strip to test the conditions.") })

  unit_info <- reactive({
    cl <- req(CELLS()); unit <- req(input$cd_unit); cond <- req(input$cd_cond)
    t <- cl$pairs[[unit]][[cond]]; req(!is.null(t)); t <- t[rowSums(t) > 0, , drop = FALSE]
    major <- setNames(colnames(t)[max.col(t, "first")], rownames(t)); purity <- apply(t, 1, max) / rowSums(t)
    keep <- rep(TRUE, length(major)); fc <- input$cd_filt_col %||% "none"
    if (fc != "none" && length(input$cd_filt_val) && !is.null(cl$pairs[[unit]][[fc]])) {
      tf <- cl$pairs[[unit]][[fc]]; fm <- setNames(colnames(tf)[max.col(tf, "first")], rownames(tf)); keep <- fm[names(major)] %in% input$cd_filt_val
    }
    list(major = major[keep], mixed = sum(purity[keep] < 0.95))
  })
  cd_base <- reactive({
    b <- req(D()); cl <- req(CELLS()); unit <- req(input$cd_unit); conds <- input$cd_conds
    validate(need(length(conds) > 0, "Tick at least one condition in the list on the left."))
    ui <- unit_info(); C <- cl$comp[[unit]]
    units <- intersect(names(ui$major)[ui$major %in% conds], colnames(C))
    validate(need(length(units) > 0, "No samples left with these choices."))
    list(C = C[, units, drop = FALSE], units = units, cond_of = factor(ui$major[units], levels = conds), conds = conds)
  })
  cd_disp <- reactive(disp_labels(req(input$cd_grp_by), cd_grp_levels()))
  cd_nfeat <- reactive(if (identical(input$cd_mode, "prop")) length(input$cd_groups) else switch(input$cd_what %||% "module", module = length(input$cd_mod), gene = length(input$cd_gene), sig = length(input$cd_sig), 1))

  cd_result <- reactive({
    b <- req(D()); base <- cd_base(); C <- base$C; units <- base$units; cond_of <- base$cond_of; conds <- base$conds
    gk <- req(input$cd_grp_by); g <- as.character(mv(gk)); validate(need(length(g) > 0, "These groups need cell information."))
    so <- style_opts("cd"); so$split <- FALSE; type <- input$cd_type %||% "box"; hc <- heat_cols()
    unit_lab <- pretty_name(input$cd_unit); cond_lab <- pretty_name(input$cd_cond); disp <- cd_disp(); test <- input$cd_test %||% "wilcox"
    if (input$cd_mode == "prop") {
      groups <- input$cd_groups; validate(need(length(groups) > 0, "Tick at least one group in the list on the left, or pick rows in the ranked table."))
      tot <- colSums(C)
      if (type == "stacked") {
        agg <- rowsum(C, ifelse(g %in% groups, g, "Other"))
        if (input$cd_stack_by == "cond") { agg <- t(rowsum(t(agg), as.character(cond_of))); xlv <- intersect(conds, colnames(agg)) } else xlv <- units[order(cond_of)]
        frac <- sweep(agg, 2, pmax(colSums(agg), 1), "/")
        glv <- c(intersect(groups, rownames(frac)), intersect("Other", rownames(frac)))
        df <- data.frame(x = factor(rep(colnames(frac), each = nrow(frac)), levels = xlv), fill = factor(rep(rownames(frac), ncol(frac)), levels = glv), y = as.vector(frac))
        df <- df[!is.na(df$fill) & !is.na(df$x), ]
        if (input$cd_stack_by == "unit") df$cond <- cond_of[match(as.character(df$x), units)]
        lab_fill <- c(disp, Other = "Other")
        p <- ggplot(df, aes(x, y, fill = fill)) + geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
          scale_fill_manual(values = c(pal_for(setdiff(glv, "Other")), Other = GREY), labels = lab_fun(lab_fill), name = key_label(gk)) +
          scale_y_continuous(labels = function(v) paste0(round(100 * v), "%"), expand = c(0, 0)) +
          labs(x = NULL, y = "Share of cells", title = paste(key_label(gk), "by", tolower(cond_lab))) + theme_mcx(so$txt) +
          theme(panel.grid.major.y = element_blank(), axis.text.x = x_text(levels(df$x), if (input$cd_stack_by == "unit") "90" else so$angle))
        if (input$cd_stack_by == "unit") p <- p + facet_grid(~ cond, scales = "free_x", space = "free_x")
        tab <- data.frame(group = unname(lab_fill[as.character(df$fill)]), column = df$x, share = round(100 * df$y, 2))
        names(tab) <- c(key_label(gk), if (input$cd_stack_by == "cond") cond_lab else unit_lab, "% of cells")
        return(list(plot = p, table = tab, stats = NULL, name = paste0(b$name, "_composition_by_", safe_name(input$cd_cond))))
      }
      rows <- lapply(groups, function(gr) data.frame(unit = units, x = cond_of, y = 100 * colSums(C[which(g == gr), , drop = FALSE]) / pmax(tot, 1), feature = gr, n = tot))
      df <- do.call(rbind, rows); df$feature <- factor(unname(disp[df$feature]), levels = unname(disp[groups])); df$fill <- df$x
      ylab <- paste0("% of the ", tolower(unit_lab), "'s cells"); stem <- "proportions"; size_by <- "mean"; size_name <- "Mean %"
      tab <- data.frame(df$unit, df$x, df$feature, round(df$y, 3), df$n); names(tab) <- c(unit_lab, cond_lab, key_label(gk), "% of cells", "Cells in sample")
      title <- NULL; sub <- NULL
    } else {
      cells <- input$cd_cells; G <- if (length(cells)) g %in% cells else rep(TRUE, length(g))
      m <- measure("cd", b$mcs); W <- C[G, , drop = FALSE]; n <- colSums(W)
      S <- as.matrix(crossprod(W, m$V[G, , drop = FALSE])) / pmax(n, 1)
      keep <- n >= (input$cd_min %||% 10)
      validate(need(any(keep), "No sample has enough cells in these groups. Lower the minimum under More options > Filters."))
      S <- S[keep, , drop = FALSE]; u <- units[keep]
      df <- data.frame(unit = rep(u, ncol(S)), x = rep(cond_of[keep], ncol(S)), feature = factor(rep(colnames(S), each = nrow(S)), levels = colnames(S)),
                       y = as.vector(S), n = rep(n[keep], ncol(S))); df$fill <- df$x
      who <- if (length(cells)) paste(unname(disp[cells]), collapse = " + ") else "all cells"
      ylab <- m$ylab; stem <- m$short; size_by <- "above"; size_name <- "% samples\nabove median"; title <- m$title; sub <- paste0("in ", who, ", averaged per ", tolower(unit_lab))
      tab <- data.frame(df$unit, df$x, df$feature, round(df$y, 4), df$n); names(tab) <- c(unit_lab, cond_lab, "Feature", ylab, "Cells used")
    }
    k <- nlevels(df$feature); nm <- paste0(b$name, "_", stem, "_by_", safe_name(input$cd_cond))
    if (type %in% c("heatmap", "dotplot")) {
      lvx <- levels(droplevels(factor(df$x))); ref <- if (isTRUE(input$cd_ref %in% lvx)) input$cd_ref else lvx[1]
      st <- if (isTRUE(input$cd_stats) && length(lvx) > 1) tile_tests(df, ref, test, input$cd_adj %||% "holm", input$cd_lab, isTRUE(input$cd_hidens)) else NULL
      p <- summary_plot(df, type, c(so, list(size_by = size_by, size_name = size_name, val_name = ylab, low = hc[1], mid = hc[2], high = hc[3],
                                             title = title, tile_labels = st, ref = if (!is.null(st)) ref, statsize = iv("cd_statsize", 3.5))))
      if (!is.null(sub)) p <- p + labs(subtitle = sub)
      return(list(plot = p, table = tab, stats = if (!is.null(st)) st[, setdiff(names(st), "label")], name = nm))
    }
    if (k > 1) df$facet <- df$feature
    p <- dist_plot(df, type, pal_for(conds), so) + labs(title = title, subtitle = sub, x = NULL, y = ylab, fill = cond_lab, colour = cond_lab)
    if (k > 1) p <- p + facet_wrap(~ facet, scales = "free_y", ncol = iv("cd_ncol", 4))
    st <- NULL
    if (isTRUE(input$cd_stats)) { o <- stat_opts("cd"); st <- run_stats(df, o); p <- add_brackets(p, st, o, if (k > 1) levels(df$feature)); st <- stats_table(st) }
    list(plot = p, table = tab, stats = st, name = nm)
  })
  output$cd_note <- renderUI({
    req(CELLS(), input$cd_unit); ui <- unit_info(); u <- tolower(pretty_name(input$cd_unit))
    txt <- if (input$cd_mode == "prop" && !identical(input$cd_type, "stacked")) paste0("Each point is one ", u, ": the share of its cells that fall in the group.")
           else if (input$cd_mode == "expr") paste0("Each point is one ", u, ": the average over its cells, where each cell takes the score of its metacell.") else NULL
    if (ui$mixed > 0) txt <- paste0(txt, " ", ui$mixed, " ", u, "(s) have cells from more than one condition and are counted under the main one.")
    if (is.null(txt)) NULL else div(class = "mcx-hint mb-1", icon("circle-info"), " ", txt)
  })
  output$cd_plot <- renderPlot({
    r <- cd_result(); reg_fig("04_conditions", "Conditions", function() print(r$plot), input[["cddl-w"]], input[["cddl-h"]])
    reg_table("conditions_values_per_sample", r$table); if (!is.null(r$stats)) reg_table("conditions_statistics", r$stats)
    r$plot
  }, height = function() heights(iv("cd_height", 480), input$cd_type %||% "box", cd_nfeat(), iv("cd_ncol", 4)), res = 96)
  dl_server("cddl", function() { p <- cd_result()$plot; function() print(p) }, function() cd_result()$name)
  output$cd_stats_tbl <- renderDT({
    validate(need(isTRUE(input$cd_stats), "Turn on Statistics above the plot."))
    s <- cd_result()$stats; validate(need(!is.null(s), "Needs at least two conditions with two or more samples.")); stats_dt(s) })
  output$cd_data_tbl <- renderDT(datatable(cd_result()$table, rownames = FALSE, selection = "none", class = "compact stripe", options = list(dom = "tp", pageLength = 10, scrollX = TRUE)))
  output$cd_dl_stats <- downloadHandler(function() paste0(cd_result()$name, "_stats.xlsx"), function(file) write_xlsx_sheets(file, list(Statistics = cd_result()$stats %||% data.frame())))
  output$cd_dl_data <- downloadHandler(function() paste0(cd_result()$name, "_values.xlsx"), function(file) write_xlsx_sheets(file, list(`Values per sample` = cd_result()$table)))
  observeEvent(input$cd_colors, {
    if (identical(input$cd_mode, "prop") && identical(input$cd_type, "stacked")) open_colors(input$cd_groups, "Group colours") else open_colors(input$cd_conds, "Condition colours")
  })

  output$cd_screen_title <- renderUI(if (identical(input$cd_mode, "expr")) "Rank all modules" else "Rank all groups")
  cd_screen <- reactive({
    b <- req(D()); base <- cd_base(); C <- base$C; cond_of <- base$cond_of; conds <- base$conds
    validate(need(length(conds) >= 2, "Tick at least two conditions to rank."))
    gk <- req(input$cd_grp_by); g <- as.character(mv(gk)); req(length(g) > 0); test <- input$cd_test %||% "wilcox"
    if (input$cd_mode == "prop") {
      tot <- pmax(colSums(C), 1); lv <- setdiff(cd_grp_levels(), "NA")
      P <- 100 * sweep(rowsum(C, g, reorder = TRUE), 2, tot, "/"); P <- P[intersect(lv, rownames(P)), , drop = FALSE]
      mean_by <- matrix(t(apply(P, 1, function(v) tapply(v, cond_of, mean))), nrow = nrow(P), dimnames = list(rownames(P), levels(cond_of)))
      pv <- apply(P, 1, function(v) p_overall(v, cond_of, test))
      d <- data.frame(key = rownames(P), Group = unname(cd_disp()[rownames(P)]), check.names = FALSE, stringsAsFactors = FALSE)
      for (j in conds) d[[paste0(j, " (mean %)")]] <- round(mean_by[, j], 2)
    } else {
      cells <- input$cd_cells; G <- if (length(cells)) g %in% cells else rep(TRUE, length(g))
      W <- C[G, , drop = FALSE]; n <- colSums(W); keep <- n >= (input$cd_min %||% 10)
      validate(need(sum(keep) >= 2, "Too few samples with enough cells in these groups."))
      S <- as.matrix(crossprod(W[, keep, drop = FALSE], b$Z[G, , drop = FALSE])) / n[keep]; cf <- droplevels(cond_of[keep])
      mean_by <- t(apply(S, 2, function(v) tapply(v, cf, mean)))
      pv <- apply(S, 2, function(v) p_overall(v, cf, test))
      i <- match(colnames(S), b$mods$id)
      d <- data.frame(key = colnames(S), Module = b$mods$short[i], Genes = vapply(b$mod_genes[colnames(S)], function(x) paste(head(x, 6), collapse = ", "), ""),
                      check.names = FALSE, stringsAsFactors = FALSE)
      for (j in intersect(conds, colnames(mean_by))) d[[paste0(j, " (mean)")]] <- round(mean_by[, j], 3)
    }
    d$`Highest in` <- colnames(mean_by)[max.col(mean_by, "first")]
    d$p <- pv; d$`Adjusted p` <- stats::p.adjust(pv, "BH")
    d[order(d$p, na.last = TRUE), , drop = FALSE]
  })
  output$cd_screen_tbl <- renderDT({
    d <- cd_screen()
    datatable(d[, -1], rownames = FALSE, selection = list(mode = "multiple"), class = "compact stripe hover nowrap",
              options = list(dom = "ftp", pageLength = 8, scrollX = TRUE, ordering = TRUE)) |> formatSignif(c("p", "Adjusted p"), 3)
  })
  observeEvent(input$cd_screen_plot, {
    d <- cd_screen(); rows <- input$cd_screen_tbl_rows_selected
    if (!length(rows)) return(warn("No rows are ticked in the ranked table.", todo = "Click rows of the table to tick them, then try again."))
    if (input$cd_mode == "prop") { if (identical(input$cd_type, "stacked")) updateRadioButtons(session, "cd_type", selected = "box"); fill_groups_rank("cd_groups", input$cd_grp_by, d$key[rows]) }
    else { updateRadioButtons(session, "cd_what", selected = "module"); set_mod_choices("cd_mod", head(d$key[rows], 20)) }
  })
  observeEvent(input$cd_screen_heat, {
    d <- cd_screen(); rows <- input$cd_screen_tbl_rows_selected
    if (!length(rows)) return(warn("No rows are ticked in the ranked table.", todo = "Click rows of the table to tick them, then try again."))
    rv$sel <- D()$mcs[mv(input$cd_grp_by) %in% d$key[rows]]; go("heat")
  })
  output$cd_dl_screen <- downloadHandler(function() paste0(D()$name, "_ranked_by_", safe_name(input$cd_cond), ".xlsx"),
                                         function(file) write_xlsx_sheets(file, list(Ranked = cd_screen()[, -1])))

  # ============================================================ HEATMAP + GENES
  g_custom_d <- debounce(reactive(input$g_custom), 700)
  gene_list <- reactive({
    b <- req(D()); src <- input$g_src %||% "top"
    switch(src,
      top = { tg <- grp_top_genes() %||% top_genes_all(b$bs, 10); unique(unlist(lapply(tg, head, input$g_top_n %||% 2))) },
      var = { v <- b$bs$gene_var[b$bs$gene_frac >= 5e-5]; head(names(sort(v, decreasing = TRUE)), input$g_var_n %||% 40) },
      mod = { id <- input$g_mod; if (isTRUE(id %in% names(b$mod_genes))) intersect(b$mod_genes[[id]], rownames(b$counts)) else character(0) },
      sig = { nm <- input$g_sig[1]; if (isTRUE(nm %in% names(rv$sigs))) rv$sigs[[nm]]$genes else character(0) },
      custom = match_genes(parse_list(g_custom_d()), rownames(b$counts))$found)
  })
  output$g_status <- renderUI({
    b <- req(D()); g <- gene_list()
    miss <- if (identical(input$g_src, "custom")) match_genes(parse_list(g_custom_d()), rownames(b$counts))$missing else character(0)
    hint(paste0(length(g), " genes"), if (length(g) > 200) " (the first 200 are drawn)", if (length(miss)) paste0(". Not found: ", paste(head(miss, 12), collapse = ", "), if (length(miss) > 12) " ..."))
  })
  add_genes <- function(g) {
    cur <- if (identical(input$g_src, "custom")) parse_list(input$g_custom) else character(0)
    updateSelectInput(session, "g_src", selected = "custom")
    updateTextAreaInput(session, "g_custom", value = paste(unique(c(cur, g)), collapse = ", "))
    go("heat2")
  }
  mc_order_file <- reactive({ f <- req(input$mc_order_file)
    x <- if (grepl("\\.xlsx$", f$name, ignore.case = TRUE)) openxlsx::read.xlsx(f$datapath, colNames = FALSE)[, 1] else utils::read.csv(f$datapath, header = FALSE, stringsAsFactors = FALSE)[, 1]
    x <- as.character(x); x[!tolower(x) %in% c("metacell", "metacells")] })
  output$mc_order_dl <- downloadHandler(function() paste0(D()$name, "_metacell_order.xlsx"), function(file) {
    gt <- grouping(); write_xlsx_sheets(file, list(Order = data.frame(metacell = gt$mc_order, block = unname(gt$id[gt$mc_order])))) })
  output$h2_title <- renderUI({ gt <- grouping()
    paste0("Columns: ", switch(gt$by, block = "pipeline blocks", cluster = paste(length(gt$levels), "new blocks"), none = "every metacell on its own",
                               ann = "names given", paste("blocks by", key_label(gt$by)))) })
  observeEvent(input$h2_go_annot, go("heat"))
  h2_cols <- reactive(c(col_or(input$h2_low, "#000080"), col_or(input$h2_mid, "#FFFFFF"), col_or(input$h2_high, "#B22222")))
  h2_data <- reactive({
    b <- req(D()); gt <- grouping(); view <- if (gt$by == "none") "mc" else input$h2_view %||% "mc"; f <- factor(gt$id, gt$levels)
    g <- head(gene_list(), 200)
    if (view == "mc") {
      cols <- gt$mc_order
      if (isTRUE(input$h2_custom_order) && !is.null(input$mc_order_file)) { cols <- intersect(mc_order_file(), cols); validate(need(length(cols) > 0, "No names in the order file match this dataset.")) }
      M <- t(b$Z[cols, , drop = FALSE]); G <- if (length(g)) as.matrix(log_expr(b, g, cols)) else NULL; grp <- unname(gt$id[cols])
    } else {
      n <- tabulate(f, nlevels(f)); cols <- gt$levels; grp <- cols
      M <- t(rowsum(b$Z, f, reorder = TRUE) / n)
      G <- if (length(g)) t(rowsum(t(as.matrix(log_expr(b, g))), f, reorder = TRUE) / n) else NULL
    }
    ro <- b$row_ord; ms <- input$h2_modsel %||% "all"
    if (ms == "pick") { keep <- intersect(input$h2_mods, b$mods$id); validate(need(length(keep) > 0, "Choose one or more modules under 'Modules on top', or show all.")); ro <- ro[b$mods$id[ro] %in% keep] }
    M <- if (ms == "none") NULL else { x <- M[ro, , drop = FALSE]; rownames(x) <- b$mods$lab[ro]; x }
    validate(need(!is.null(M) || !is.null(G), "Nothing to draw: no modules on top and no genes underneath."))
    if (!is.null(G)) { rownames(G) <- g; if (isTRUE(input$h2_genez)) { G <- t(scale(t(G))); G[is.na(G)] <- 0 } }
    bb <- label_bars(input$h2_bars, view, cols, f)
    list(M = M, G = G, grp = if (gt$by == "none") rep("x", length(grp)) else grp, bars = bb$bars, bc = bb$bc)
  })
  output$h2_legend <- renderUI(legend_ui(h2_data()$bc))
  observeEvent(input$h2_colors, { d <- h2_data()
    open_colors(unlist(lapply(d$bc, names)), "Heatmap colours", scale = list(ids = c("h2_low", "h2_mid", "h2_high"), def = c("#000080", "#FFFFFF", "#B22222"),
      extra = tagList(div(class = "mcx-grid2 mt-1",
        num_step("h2_lim", "Module colour range (z)", iv("h2_lim", 2), 0.5, 5, 0.25),
        selectInput("h2_gpal", lab("Gene colours", "Colours of the gene rows, from no expression to the highest."),
                    c("Viridis" = "viridis", "White to red" = "reds", "White to blue" = "blues", "Magma" = "magma", "Same as modules" = "same"), selected = input$h2_gpal %||% "viridis"))))) })
  h2_plot <- reactive({
    d <- h2_data(); r <- rle(as.character(d$grp)); ends <- cumsum(r$lengths); mids <- ends - r$lengths / 2 + 0.5
    wide <- r$lengths >= max(1, length(d$grp) * 0.012)
    top <- if (isTRUE(input$h2_sep) && sum(wide) <= 60 && length(r$values) > 1) list(at = mids[wide], lab = r$values[wide]) else NULL
    ht <- iv("h2_height", 700); nrow_all <- (if (is.null(d$M)) 0 else nrow(d$M)) + (if (is.null(d$G)) 0 else nrow(d$G))
    ysize <- max(4, min(9, 0.6 * ht / max(nrow_all, 1)))
    hc <- h2_cols(); lim <- iv("h2_lim", 2)
    parts <- list(); hts <- numeric(0)
    if (length(d$bars)) { parts <- c(parts, list(bar_panel(d$bars, d$bc, ysize))); hts <- c(hts, 0.03 * length(d$bars)) }
    gp <- if (is.null(d$G)) 0 else if (is.null(d$M)) 1 else iv("h2_gprop", 0.6)
    if (gp > 0) {   # genes first, right under the label bars
      pal <- input$h2_gpal %||% "viridis"
      gcol <- if (isTRUE(input$h2_genez) || pal == "same") hc else switch(pal, reds = c("#FFFFFF", "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D"),
                blues = c("#FFFFFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"), magma = grDevices::hcl.colors(9, "Inferno"), grDevices::hcl.colors(9, "viridis"))
      glim <- if (isTRUE(input$h2_genez)) c(-lim, lim) else c(0, max(stats::quantile(d$G, 0.99), 1e-6))
      parts <- c(parts, list(heat_panel(d$G, gcol, glim, if (isTRUE(input$h2_genez)) "Gene (z)" else "Expression\nlog2(1+UMI/100k)", ysize, ends, isTRUE(input$h2_sep),
                                        if (!length(d$bars)) top)))
      hts <- c(hts, gp)
    }
    if (!is.null(d$M)) { parts <- c(parts, list(heat_panel(d$M, hc, c(-lim, lim), "Module\nscore (z)", ysize, ends, isTRUE(input$h2_sep), if (gp == 0 && !length(d$bars)) top))); hts <- c(hts, 1 - gp) }
    wrap_plots(parts, ncol = 1, heights = hts)
  })
  output$h2_plot <- renderPlot({ p <- h2_plot(); reg_fig("02_heatmap_genes", "Heatmap & genes", function() print(p), input[["h2dl-w"]], input[["h2dl-h"]]); p },
                               height = function() iv("h2_height", 760), res = 96)
  dl_server("h2dl", function() { p <- h2_plot(); function() print(p) }, function() paste0(D()$name, "_heatmap_genes"))

  # ============================================================ MODULES & MARKERS
  by_group_plot <- function(v, key, title, ylab) {
    x <- as.character(mv(key)); validate(need(length(x) > 0, "These groups need cell information.")); x[is.na(x)] <- "NA"
    lv <- level_order(key, x); validate(need(length(lv) <= 400, "Too many groups to plot; choose another grouping."))
    df <- data.frame(x = factor(x, levels = lv), y = v); df$fill <- df$x
    dist_plot(df, "box", pal_for(lv), list(points = TRUE, pt = 0.4, alpha = 0.3, txt = 12, legend = FALSE, angle = if (length(lv) > 12) "90" else "auto",
                                          xlabels = disp_labels(key, lv))) + labs(title = title, x = NULL, y = ylab)
  }
  mod_table <- reactive({
    b <- req(D()); key <- req(input$mg_x); x <- as.character(mv(key)); req(length(x) > 0); x[is.na(x)] <- "NA"
    lv <- level_order(key, x); f <- factor(x, lv); o <- order(b$mods$num)
    zm <- rowsum(b$Z, f, reorder = TRUE) / tabulate(f, length(lv)); dl <- disp_labels(key, lv)
    top3 <- apply(zm, 2, function(z) { k <- order(z, decreasing = TRUE)[1:min(3, length(z))]; paste0(dl[rownames(zm)[k]], " (", round(z[k], 1), ")", collapse = "; ") })
    data.frame(id = b$mods$id[o], Module = b$mods$short[o], Genes = b$mods$n_genes[o], `Highest in` = top3[b$mods$id[o]], `Gene list` = b$mods$genes[o], check.names = FALSE, stringsAsFactors = FALSE)
  })
  output$mod_tbl <- renderDT({
    d <- mod_table()
    datatable(d[, -1], rownames = FALSE, selection = "single", class = "compact stripe hover",
              options = list(dom = "ftp", pageLength = 8, scrollX = TRUE, language = list(search = "Find a module or gene:"),
                             columnDefs = list(list(className = "wrap", targets = 3))))
  })
  observeEvent(input$mod_tbl_rows_selected, { d <- mod_table(); r <- input$mod_tbl_rows_selected; if (length(r)) set_mod_choices("fm_mod", d$id[r]) })
  fm_id <- reactive({ b <- req(D()); id <- input$fm_mod; if (!isTRUE(id %in% b$mods$id)) id <- b$mods$id[order(b$mods$num)][1]; id })
  output$fm_title <- renderUI({ b <- req(D()); i <- match(fm_id(), b$mods$id); paste0(b$mods$name[i], "  ·  ", b$mods$n_genes[i], " genes") })
  output$fm_genes <- renderUI(genes_ui(D()$mod_genes[[fm_id()]], 80, copy = TRUE))
  fm_plot <- reactive({ b <- req(D()); id <- fm_id(); i <- match(id, b$mods$id); by_group_plot(b$Z[, id], req(input$mg_x), b$mods$title[i], "Module score (z)") })
  output$fm_plot <- renderPlot({ p <- fm_plot(); reg_fig("05_module", "Modules & markers: module", function() print(p), input[["fmdl-w"]], input[["fmdl-h"]]); p }, res = 96)
  dl_server("fmdl", function() { p <- fm_plot(); function() print(p) }, function() paste0(D()$name, "_", D()$mods$short[match(fm_id(), D()$mods$id)]))
  output$mod_dl <- downloadHandler(function() paste0(D()$name, "_gene_modules.xlsx"), function(file) write_xlsx_sheets(file, list(Modules = mod_table()[, -1])))
  feature_to <- function(to, what, items, key = NULL) {
    updateRadioButtons(session, paste0(to, "_what"), selected = what)
    if (what == "module") set_mod_choices(paste0(to, "_mod"), items) else if (what == "gene") set_gene_choices(paste0(to, "_gene"), items)
    if (to == "pl" && !is.null(key)) updateSelectInput(session, "pl_x", selected = key)
  }
  observeEvent(input$fm_go_plot, { feature_to("pl", "module", fm_id(), input$mg_x); go("plot") })
  observeEvent(input$fm_go_cond, { if (is.null(CELLS())) return(warn("Conditions need cell information.", todo = "Add the Seurat object or a table of cells on the Home page.")); feature_to("cd", "module", fm_id()); updateRadioButtons(session, "cd_mode", selected = "expr"); go("cond") })
  observeEvent(input$fm_sig, { b <- D(); open_sig_modal(NULL, b$mod_genes[[fm_id()]], b$mods$name[match(fm_id(), b$mods$id)]) })
  observeEvent(input$fm_go_h2, { updateSelectInput(session, "g_src", selected = "mod"); set_mod_choices("g_mod", fm_id()); go("heat2") })
  observeEvent(input$fm_go_umap, { updateSelectInput(session, "um_by", selected = "module"); set_mod_choices("um_mod", fm_id()); go("umap") })

  output$fg_title <- renderUI(req(input$fg_gene))
  output$fg_info <- renderUI({
    b <- req(D()); g <- req(input$fg_gene); m <- b$gene2mod[g]
    if (is.na(m)) return(div(class = "mcx-hint", strong(g), " is not in any gene module (not variable enough to be used for modules)."))
    i <- match(m, b$mods$short)
    tagList(div(class = "mcx-hint mb-1", strong(g), " is in ", strong(b$mods$name[i]), paste0(" (", b$mods$n_genes[i], " genes):")), genes_ui(b$mod_genes[[b$mods$id[i]]], 40))
  })
  fg_plot <- reactive({ b <- req(D()); g <- req(input$fg_gene); req(g %in% rownames(b$counts)); by_group_plot(gene_expr(b, g), req(input$mg_x), g, "Expression, log2(1 + UMIs per 100k)") })
  output$fg_plot <- renderPlot({ p <- fg_plot(); reg_fig("06_gene", "Modules & markers: gene lookup", function() print(p), input[["fgdl-w"]], input[["fgdl-h"]]); p }, res = 96)
  dl_server("fgdl", function() { p <- fg_plot(); function() print(p) }, function() paste0(D()$name, "_", input$fg_gene))
  observeEvent(input$fg_go_plot, { feature_to("pl", "gene", input$fg_gene, input$mg_x); go("plot") })
  observeEvent(input$fg_go_cond, { if (is.null(CELLS())) return(warn("Conditions need cell information.", todo = "Add the Seurat object or a table of cells on the Home page.")); feature_to("cd", "gene", input$fg_gene); updateRadioButtons(session, "cd_mode", selected = "expr"); go("cond") })
  observeEvent(input$fg_go_umap, { updateSelectInput(session, "um_by", selected = "gene"); set_gene_choices("um_gene", input$fg_gene); go("umap") })
  observeEvent(input$fg_go_mod, { b <- D(); m <- b$gene2mod[input$fg_gene]; if (is.na(m)) return(warn(input$fg_gene, " is not in any module.", todo = "Its expression is still shown in the plot."))
    set_mod_choices("fm_mod", b$mods$id[match(m, b$mods$short)]); updateRadioButtons(session, "mg_view", selected = "mods") })

  tg_levels <- reactive({ req(D(), input$tg_from); group_levels(input$tg_from) })
  observeEvent(tg_levels(), {
    lv <- tg_levels(); dl <- disp_labels(input$tg_from, lv); big <- length(lv) > 300
    pend <- isolate(rv$pending$tg_a); rv$pending$tg_a <- NULL
    updateSelectizeInput(session, "tg_a", choices = setNames(lv, dl[lv]), selected = pend %||% lv[1], server = big)
    updateSelectizeInput(session, "tg_b", choices = setNames(lv, dl[lv]), selected = character(0), server = big)
  })
  output$tg_title <- renderUI({ a <- input$tg_a; if (!length(a)) "Marker genes" else paste("Marker genes of", paste(unname(disp_labels(input$tg_from, a)), collapse = " + "),
                                                                                          if (identical(input$tg_vs, "groups") && length(input$tg_b)) paste("vs", paste(input$tg_b, collapse = " + ")) else "vs all other metacells") })
  tg_result <- reactive({
    b <- req(D()); key <- req(input$tg_from); A <- input$tg_a; validate(need(length(A) > 0, "Choose a group on the left."))
    if (key == "grp") { gt <- grouping(); validate(need(length(gt$levels) <= 400, "Too many heatmap groups; choose a coarser grouping on Annotate.")); gs <- grouping_gstats(); glv <- gt$levels }
    else {
      v <- as.character(mv(key)); v[is.na(v)] <- "NA"; glv <- level_order(key, v)
      validate(need(length(glv) <= 400, "Too many groups."))
      gs <- if (key == "block") b$bs else gstats_for(v, glv, paste(key, length(glv), sum(match(v, glv) * (seq_along(v) %% 97))))
    }
    B <- if (identical(input$tg_vs, "groups")) setdiff(input$tg_b %||% character(0), A) else NULL
    validate(need(!identical(input$tg_vs, "groups") || length(B), "Choose groups to compare with."))
    r <- head(rank_genes(gs, match(A, glv), if (length(B)) match(B, glv) else NULL), input$tg_n %||% 50)
    data.frame(Gene = r$gene, `log2 fold change` = round(r$log2FC, 2), `Mean fraction` = signif(r$mean_frac, 3), Module = unname(b$gene2mod[r$gene]), check.names = FALSE)
  })
  output$tg_tbl <- renderDT({ d <- tg_result(); reg_table("marker_genes", d)
    datatable(d, rownames = FALSE, selection = "none", class = "compact stripe", options = list(pageLength = 15, dom = "ftp"),
              container = dt_head(names(d), c(Gene = "Gene symbol", `log2 fold change` = "Mean expression in the group against the comparison, log2. Positive: higher in the group.",
                `Mean fraction` = "Average share of the group's counts that come from this gene.", Module = "Gene module the gene belongs to, if any."))) })
  output$tg_dl <- downloadHandler(function() paste0(D()$name, "_marker_genes_", safe_name(paste(input$tg_a, collapse = "_")), ".xlsx"), function(file) write_xlsx_sheets(file, list(`Marker genes` = tg_result())))
  observeEvent(input$tg_go_h2, add_genes(head(tg_result()$Gene, 20)))
  observeEvent(input$tg_sig, open_sig_modal(NULL, head(tg_result()$Gene, 20), paste("Markers", paste(input$tg_a, collapse = "+"))))
  observeEvent(input$tg_go_plot, { feature_to("pl", "gene", head(tg_result()$Gene, 6), input$tg_from)
    if (identical(input$pl_x, input$tg_from)) fill_groups_rank("pl_groups", input$tg_from) ; go("plot") })

  vg_levels <- reactive({ req(D(), input$vg_from); group_levels(input$vg_from) })
  observeEvent(vg_levels(), { lv <- vg_levels(); lv <- if (length(lv) > 400) character(0) else lv
    want <- rv$pending$vg_groups; rv$pending$vg_groups <- NULL
    updateSelectizeInput(session, "vg_groups", choices = setNames(lv, disp_labels(input$vg_from, lv)[lv]), selected = intersect(want %||% isolate(input$vg_groups), lv)) })
  output$vg_note <- renderUI({ if (is.null(CELLS()$total_umis)) hint(icon("circle-info"), " Needs cell information with UMI counts (nCount_RNA or total_umis).") })
  observeEvent(input$vg_go, {
    b <- req(D()); cl <- CELLS(); g <- input$vg_groups
    if (!length(g)) return(warn("No groups are chosen.", todo = "Pick one or more groups in the box above Find genes."))
    if (is.null(cl$total_umis)) return(warn("This needs UMI counts per cell (nCount_RNA) in the cell information.", todo = "Use a Seurat object with nCount_RNA."))
    mcs <- b$mcs[mv(input$vg_from) %in% g]; if (length(mcs) < 3) return(warn("The chosen groups have fewer than 3 metacells.", todo = "Add more groups."))
    withProgress(message = "Testing genes", {
      umi <- cl$total_umis[match(mcs, b$mcs)]
      O <- as(b$counts[, mcs, drop = FALSE] %*% Diagonal(x = umi / b$libsize[mcs]), "CsparseMatrix")
      rs <- Matrix::rowSums(O); mx <- numeric(nrow(O))
      if (length(O@x)) { v <- tapply(O@x, O@i + 1L, max); mx[as.integer(names(v))] <- v }
      keep <- mx > 3 & rs >= 10; O <- O[keep, , drop = FALSE]; rs <- rs[keep]
      p <- umi / sum(umi); O2 <- O; O2@x <- O2@x^2
      chi <- as.vector(O2 %*% (1 / p)) / rs - rs
      pv <- stats::pchisq(chi, df = length(mcs) - 1, lower.tail = FALSE)
      res <- data.frame(Gene = rownames(O), `Chi-square` = round(chi, 1), p = signif(pv, 3), `Adjusted p` = signif(stats::p.adjust(pv, "BH"), 3),
                        UMIs = round(rs), Module = unname(b$gene2mod[rownames(O)]), check.names = FALSE, stringsAsFactors = FALSE)
      rv$vg <- res[order(-res$`Chi-square`), ]; reg_table("genes_varying_within_groups", rv$vg)
    })
  })
  output$has_vg <- reactive(!is.null(rv$vg)); outputOptions(output, "has_vg", suspendWhenHidden = FALSE)
  output$vg_tbl <- renderDT({ validate(need(!is.null(rv$vg), "Choose one or more blocks on the left and click Find varying genes."))
    d <- head(rv$vg, 1000)
    datatable(d, rownames = FALSE, selection = "none", class = "compact stripe", options = list(pageLength = 15, dom = "ftp"),
              container = dt_head(names(d), c(Gene = "Gene symbol", `Chi-square` = "How much the gene's counts differ between the block's metacells beyond what their sizes explain. Higher means more variation.",
                p = "Chance of this much variation if the gene were equally expressed in every metacell of the block.", `Adjusted p` = "p corrected for testing every gene (Benjamini-Hochberg).",
                UMIs = "Total counts of the gene in the chosen blocks; low counts make the test less reliable.", Module = "Gene module the gene belongs to, if any."))) })
  output$vg_dl <- downloadHandler(function() paste0(D()$name, "_genes_varying_within_groups.xlsx"), function(file) write_xlsx_sheets(file, list(`Varying genes` = req(rv$vg))))
  observeEvent(input$vg_go_h2, { req(rv$vg); add_genes(head(rv$vg$Gene[rv$vg$`Adjusted p` < 0.01], 30)) })
  observeEvent(input$vg_go_annot, { g <- input$vg_groups; if (length(g) && identical(input$vg_from, "grp")) rv$sel <- group_mcs(g); go("heat") })

  # ============================================================ UMAP
  compute_umap <- function() {
    b <- req(D())
    if (!requireNamespace("uwot", quietly = TRUE)) { rv$umap_msg <- "The UMAP needs the 'uwot' package: install.packages(\"uwot\")"; return() }
    withProgress(message = "Making a UMAP from the module scores", value = 0.3, {
      set.seed(42)
      xy <- uwot::umap(b$Z, n_neighbors = min(input$um_nn, nrow(b$Z) - 1), min_dist = input$um_md, n_threads = max(1, min(4, parallel::detectCores() - 1)))
      u <- list(xy = data.frame(x = xy[, 1], y = xy[, 2], row.names = b$mcs),
                source = sprintf("Made from the module scores (%d neighbours, minimum distance %.2f)", input$um_nn, input$um_md))
      UMAP(u); cache_write(b$cache_file, b$cache_key, "umap", u); rv$umap_msg <- NULL
    })
  }
  observeEvent(input$nav, { if (input$nav == "umap" && !is.null(D()) && is.null(UMAP())) compute_umap() })
  observeEvent(input$um_recompute, compute_umap())
  output$um_source <- renderUI(hint(UMAP()$source %||% rv$umap_msg %||% ""))
  um_feature <- reactive({
    b <- req(D()); by <- req(input$um_by)
    if (by == "module") { id <- req(input$um_mod); req(id %in% b$mods$id); return(list(v = b$Z[, id], title = b$mods$title[match(id, b$mods$id)], num = TRUE, div = TRUE)) }
    if (by == "gene") { g <- req(input$um_gene); req(g %in% rownames(b$counts)); return(list(v = gene_expr(b, g), title = g, num = TRUE, div = FALSE)) }
    if (by == "sig") { nm <- input$um_sig; validate(need(!is.null(rv$sigs[[nm %||% ""]]), "Make a signature first: click New on the left.")); return(list(v = sig_scores(nm), title = nm, num = TRUE, div = TRUE)) }
    v <- mv(by); validate(need(!is.null(v), "This needs cell information.")); list(v = v, title = key_label(by), num = is.numeric(v), div = FALSE, key = by)
  })
  output$um_title <- renderUI(um_feature()$title)
  um_plot <- reactive({
    b <- req(D()); u <- UMAP(); validate(need(!is.null(u), rv$umap_msg %||% "Making the UMAP..."))
    f <- um_feature(); xy <- u$xy; xy$v <- unname(f$v[match(rownames(xy), b$mcs)])
    if (f$num) xy <- xy[order(xy$v), ] else { lv <- level_order(f$key, xy$v); xy$v <- factor(xy$v, levels = lv) }
    p <- ggplot(xy, aes(x, y, colour = v)) + geom_point(size = iv("um_pt", 1.2), stroke = 0) + coord_equal(clip = "off") +
      scale_x_continuous(expand = expansion(mult = 0.06)) + scale_y_continuous(expand = expansion(mult = 0.06)) +
      labs(x = "UMAP 1", y = "UMAP 2", colour = NULL) + theme_mcx() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.major.y = element_blank(),
            axis.line.x = element_blank(), panel.border = element_rect(fill = NA, colour = "#E5E7EB"), plot.margin = margin(10, 14, 6, 6))
    if (f$num) {
      pal <- input$um_pal %||% "div"; hc <- heat_cols()
      p <- p + if (pal == "div" && f$div) scale_colour_gradient2(low = hc[1], mid = "#E9EAEE", high = hc[3], limits = c(-2, 2), oob = scales::squish)
               else if (pal == "magma") scale_colour_viridis_c(option = "magma") else scale_colour_viridis_c()
    } else {
      p <- p + scale_colour_manual(values = pal_for(lv), na.value = GREY, labels = lab_fun(disp_labels(f$key, lv))) + guides(colour = guide_legend(override.aes = list(size = 3.5)))
      if (isTRUE(input$um_labels) && length(lv) <= 60) {
        cen <- stats::aggregate(cbind(x, y) ~ v, data = xy, FUN = stats::median)
        no_border <- if (utils::packageVersion("ggplot2") >= "3.5.0") list(linewidth = 0) else list(label.size = 0)   # label.size was renamed in ggplot2 3.5
        p <- p + do.call(geom_label, c(list(data = cen, mapping = aes(x, y, label = v), inherit.aes = FALSE, size = 3, fill = scales::alpha("white", 0.8), colour = "#1D2129"), no_border))
      }
    }
    p
  })
  output$um_plot <- renderPlot({ p <- um_plot(); reg_fig("07_umap", "UMAP", function() print(p), input[["umdl-w"]], input[["umdl-h"]]); p },
                               height = function() iv("um_height", 560), res = 96)
  dl_server("umdl", function() { p <- um_plot(); function() print(p) }, function() paste0(D()$name, "_umap_", safe_name(input$um_by)))
  observeEvent(input$um_colors, { f <- um_feature(); if (f$num) return(warn("Scores use a colour scale, not group colours.", todo = "Change it under Style on the left."))
    open_colors(level_order(f$key, f$v)) })

  # ============================================================ EXPORT
  observe({ s <- rv$figs_seen; lab <- vapply(s, function(k) get(k, envir = figs)$label, ""); updateCheckboxGroupInput(session, "ex_figs", choices = setNames(s, lab), selected = s) })
  observe({ s <- rv$tables_seen; updateCheckboxGroupInput(session, "ex_tables", choices = setNames(s, sub("^(.)", "\\U\\1", pretty_name(s), perl = TRUE)), selected = s) })
  output$ex_figs_none <- renderUI(if (!length(rv$figs_seen)) hint("No figures opened yet."))
  session_obj <- function() list(app = "MCexplorer", version = APP_VERSION, saved = Sys.time(), dataset = D()$name,
    ann = rv$ann, log = rv$log, colors = rv$colors, gops = rv$gops, split_sig = rv$split_sig, sigs = rv$sigs, pairs = rv$pairs,
    inputs = lapply(setNames(names(RESTORE), names(RESTORE)), function(id) isolate(input[[id]])))
  apply_session <- function(s) {
    rv$sigs <- s$sigs %||% list(); rv$pairs <- s$pairs %||% list(pl = character(0), cd = character(0))
    for (id in names(s$inputs)) {
      v <- s$inputs[[id]]; type <- RESTORE[[id]]
      if (is.null(type) || (is.null(v) && type != "multi")) next
      switch(type, text = updateTextAreaInput(session, id, value = v), numeric = updateNumericInput(session, id, value = v),
             slider = updateSliderInput(session, id, value = v), switch = update_switch(id, value = v),
             radio = updateRadioButtons(session, id, selected = v), select = updateSelectInput(session, id, selected = v),
             colour = colourpicker::updateColourInput(session, id, value = v),
             multi = updateSelectizeInput(session, id, selected = v %||% character(0)))
    }
    a <- rv$ann; i <- match(a$metacell, s$ann$metacell); a$annotation[!is.na(i)] <- s$ann$annotation[i[!is.na(i)]]
    rv$ann <- a; rv$log <- rbind(data.frame(time = format(Sys.time(), "%H:%M:%S"), action = "session restored", metacells = sum(!is.na(i)), name = ""), s$log %||% NULL)
    rv$colors <- s$colors %||% character(0); rv$gops <- s$gops %||% lapply(names(s$splits), function(G) list(type = "split", g = G, n = s$splits[[G]])); rv$split_sig <- s$split_sig %||% ""
  }
  observeEvent(input$session_file, {
    s <- tryCatch(readRDS(input$session_file$datapath), error = function(e) NULL)
    if (is.null(s) || !identical(s$app, "MCexplorer")) { rv$session_msg <- "This is not an MCexplorer session file."; return() }
    if (!identical(s$dataset, D()$name)) { rv$session_msg <- paste0("This session belongs to '", s$dataset, "'. That dataset has to be open first."); return() }
    apply_session(s); rv$session_msg <- paste0("Restored the session saved on ", format(s$saved, "%d %b %Y %H:%M"), ".")
  })
  output$session_msg <- renderUI(if (!is.null(rv$session_msg)) div(class = "alert alert-info py-2 px-3 small mt-2", rv$session_msg))
  output$dl_session <- downloadHandler(function() paste0(D()$name, "_session_", format(Sys.time(), "%Y%m%d_%H%M"), ".rds"), function(file) saveRDS(session_obj(), file))
  output$dl_all <- downloadHandler(
    function() paste0(D()$name, "_MCexplorer_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip"),
    function(file) {
      td <- tempfile("mcx_"); dir.create(td); fmt <- input$ex_fmt %||% "pdf"
      withProgress(message = "Preparing the zip", {
        if (isTRUE(input$ex_ann)) write_xlsx_sheets(file.path(td, "annotation.xlsx"), ann_sheets())
        if (length(input$ex_tables)) { tb <- list(); for (n in input$ex_tables) if (exists(n, envir = tables)) tb[[make.unique(c(names(tb), substr(n, 1, 31)))[length(tb) + 1]]] <- get(n, envir = tables)
          if (length(tb)) write_xlsx_sheets(file.path(td, "tables.xlsx"), tb) }
        incProgress(0.3, detail = "figures")
        if (length(input$ex_figs)) { dir.create(file.path(td, "figures"))
          for (n in input$ex_figs) { f <- get(n, envir = figs)
            tryCatch(save_plot(file.path(td, "figures", paste0(n, ".", fmt)), f$draw, f$w, f$h, fmt), error = function(e) mcx_log("Figure ", n, " failed: ", conditionMessage(e))) } }
        files <- list.files(td, recursive = TRUE)
        if (!length(files)) { writeLines("Nothing was selected.", file.path(td, "empty.txt")); files <- "empty.txt" }
        if (requireNamespace("zip", quietly = TRUE)) zip::zip(file, files, root = td)
        else { owd <- setwd(td); on.exit(setwd(owd)); utils::zip(file, files) }
      })
    })

  session$onSessionEnded(function() { j <- isolate(rv$job); if (!is.null(j$job) && j$job$is_alive()) j$job$kill() })
}
