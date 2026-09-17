# Sourced at the top of every Rmd. Nothing below this line is analysis, it is
# only the paths and factor levels that all the steps have to agree on.
#
# Usage from a script in scripts/:  source(here::here("config", "config.R"))

# repo root. Walks up from the working directory until it finds config/, so the
# Rmds knit correctly whether they are run from scripts/ or from the repo root.
if (!exists("PROJ")) {
  PROJ <- normalizePath(".", mustWork = TRUE)
  while (!dir.exists(file.path(PROJ, "config")) && dirname(PROJ) != PROJ)
    PROJ <- dirname(PROJ)
  stopifnot(dir.exists(file.path(PROJ, "config")))
}

PROJECT_NAME <- "project"

# output trees. Each step writes into its own <slug>/metacell subfolder so the
# global run and every subset run can coexist without overwriting each other.
DIR_H5AD <- file.path(PROJ, "results", "h5ad")
DIR_RDS  <- file.path(PROJ, "results", "rds")
DIR_CSV  <- file.path(PROJ, "results", "csvs")
DIR_PDF  <- file.path(PROJ, "results", "pdfs")

# sample sheet: one row per sample, the only place the design is written down
SAMPLES <- read.csv(file.path(PROJ, "config", "samples.csv"), stringsAsFactors = FALSE)

# factor levels. Declared here rather than inferred from the data so that a
# missing timepoint in one subset cannot silently reorder an axis.
COND_LEVELS <- c("CTRL", "CASE")
TP_LEVELS   <- c("T0", "T1", "T2", "T3", "T4")

# not cell types: MC2 outliers, metacells that never reached the module heatmap,
# and blocks whose lineage score was too weak to call
drop_labels <- c("Outlier", "Unassigned", "Undefined")

# the metadata column holding the final per-cell annotation
LABEL_COL <- "cluster_label"

# shared palettes, keyed by name so adding a level cannot shift existing colours
PAL_COND <- c(CTRL = "#4C7FB8", CASE = "#A02C2C")
PAL_TP   <- c(T0 = "#D9D9D9", T1 = "#A6BDDB", T2 = "#67A9CF",
              T3 = "#3690C0", T4 = "#02648A")

# helper: build an output directory tree for one analysis slug and return it
out_dirs <- function(slug, kind) {
  d <- file.path(switch(kind, csv = DIR_CSV, pdf = DIR_PDF, rds = DIR_RDS,
                        h5ad = DIR_H5AD), slug, "metacell")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
