"""
01 Preprocess: load per-sample matrices, collapse duplicated gene symbols,
remove TCR-defined doublets, attach design metadata, merge, write QC.

Output: one merged .h5ad plus a QC workbook. Nothing here is metacell-specific;
step 02 picks up the merged object.

Background and the reasoning behind each filter: docs/03_preprocessing.md
"""

import sys
print("PYTHON EXECUTABLE:", sys.executable)
print("PYTHON VERSION:", sys.version)

# ---------------------------------------------------------------------------
# 1. Libraries
# ---------------------------------------------------------------------------
import scanpy as sc            # 10x h5 reader
import anndata as ad           # AnnData containers and concat
import metacells as mc         # imported here so a broken env fails fast, used from 02
import numpy as np
import pandas as pd
import os
import glob
import yaml

# ---------------------------------------------------------------------------
# 2. Configuration
# ---------------------------------------------------------------------------
# Everything comes from config/, so this script has no hardcoded sample IDs
# and can be pointed at a different study by editing two csv files.
CFG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config")
CFG = yaml.safe_load(open(os.path.join(CFG_DIR, "config.yaml")))
ROOT = os.path.abspath(os.path.join(CFG_DIR, ".."))

def rp(p):
    """Resolve a config path against the repo root unless it is already absolute."""
    return p if os.path.isabs(p) else os.path.join(ROOT, p)

PROJECT   = CFG["project"]["name"]
data_dir  = rp(CFG["paths"]["gex_h5"])
h5ad_dir  = rp(os.path.join(CFG["paths"]["h5ad"], "global", "metacell"))
qc_csv_dir = rp(os.path.join(CFG["paths"]["csvs"], "global", "qc", "metacell"))
os.makedirs(h5ad_dir, exist_ok=True)
os.makedirs(qc_csv_dir, exist_ok=True)

# Sample sheet: one row per sample, carrying the whole design. Reading it here
# means the metadata block further down is a join, not a set of literal dicts.
samples_df = pd.read_csv(rp(CFG["paths"]["samples"]), dtype=str)
samples_df["timepoint_days"] = samples_df["timepoint_days"].astype(int)
sampleNames = samples_df["sample_origin"].tolist()
print(f"{len(sampleNames)} samples to load")

# Ordered factor levels, so every downstream plot and model sees the same order
condition_levels = CFG["design"]["condition_levels"]
category_levels = sorted(samples_df["category"].unique())
category_plot_levels = samples_df.sort_values(
    ["category", "timepoint"])["category_plot"].tolist()

# V(D)J configuration, consumed by step 5c below
VDJ_RAW_DIR = rp(CFG["paths"]["vdj_root"])
VDJ_BATCHES = CFG["paths"]["vdj_batches"]

# library folder token -> sample id. Library naming is rarely derivable from the
# sample id by pattern, so the mapping is written out explicitly in a csv.
vdj_map_df = pd.read_csv(rp(CFG["paths"]["vdj_map"]), dtype=str)
VDJ_TOKEN_TO_SAMPLE = dict(zip(vdj_map_df["library_token"], vdj_map_df["sample_origin"]))
assert set(VDJ_TOKEN_TO_SAMPLE.values()) <= set(sampleNames), \
    "vdj_library_map.csv references a sample that is not in samples.csv"

# a second beta must be independently supported, not an ambient fragment
MIN_SECOND_BETA_UMIS = CFG["preprocess"]["min_second_beta_umis"]
MIN_SECOND_BETA_FRAC = CFG["preprocess"]["min_second_beta_frac"]

# QC accumulators: one row per sample
qc_rows = []                 # summary metrics per sample
umi_pct_rows = []            # UMI distribution per sample (percentile table)
PERCENTILES = [0, 1, 5, 10, 25, 50, 75, 90, 95, 99]

#------------------------------------------

#3# Load each sample's filtered feature-barcode matrix
print("Loading datasets...")
adatas = []
for sample in sampleNames:
    filtered_file = f"{data_dir}/{sample}_filtered_feature_bc_matrix.h5"
    print(f"  Loading {sample} from {filtered_file}")
    adatas.append(sc.read_10x_h5(filtered_file))

print("Datasets loaded successfully")

#------------------------------------------

#4#data distribution per sample
print("Exploring data distribution...")
for i, adata in enumerate(adatas):

    summed_expression = np.array(adata.X.sum(axis=1)).flatten()

    print(f"\nSample: {sampleNames[i]}")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")

    print("Min UMI:", summed_expression.min())
    print("Median UMI:", np.median(summed_expression))
    print("Max UMI:", summed_expression.max())

    values = np.percentile(summed_expression, PERCENTILES)
    print("UMI percentiles:")
    for p, v in zip(PERCENTILES, values):
        print(f"  {p:>2}% : {v:.1f}")

    #append QC file
    row = {
        "sample_origin": sampleNames[i],
        "n_cells_loaded": int(adata.n_obs),
        "n_genes_loaded": int(adata.n_vars),
        "umi_min": float(summed_expression.min()),
        "umi_median": float(np.median(summed_expression)),
        "umi_max": float(summed_expression.max()),
        "umi_mean": float(summed_expression.mean()),
    }
    qc_rows.append(row)
    umi_pct_rows.append(
        {"sample_origin": sampleNames[i], **{f"p{p}": float(v) for p, v in zip(PERCENTILES, values)}}
    )

# index qc_rows by sample so later steps can update the same dict
qc_by_sample = {r["sample_origin"]: r for r in qc_rows}

#---------------------------------------

#5# Identify duplicated genes and merge them
for i in range(len(adatas)):
    adata = adatas[i].copy()
    print("Before:", adata.shape)
    gene_names = adata.var_names
    duplicated_genes = gene_names[gene_names.duplicated()].unique()
    print(f"Sample {sampleNames[i]}: {len(duplicated_genes)} duplicated genes")

    for gene in duplicated_genes:
        # identifies gene column indices
        gene_mask = (adata.var_names == gene)
        # sums expression values cell by cell
        summed_expression = np.array(adata[:, gene_mask].X.sum(axis=1)).flatten()
        # takes metadata from the gene
        gene_var = adata[:, gene_mask].var.iloc[0]
        # erases all the copies
        adata = adata[:, ~gene_mask]
        # replaces by one column with pooled values
        new_col = ad.AnnData(X=summed_expression[:, None], obs=adata.obs.copy(), var=pd.DataFrame(index=[gene], data=gene_var).T)

        adata = ad.concat([adata, new_col], axis=1)
    adatas[i] = adata
    print("After:", adata.shape)

# Check that there are no duplicates left:
for i in range(len(adatas)):
    adata = adatas[i]
    gene_names = adata.var_names
    duplicated_genes = gene_names[gene_names.duplicated()].unique()
    print(sampleNames[i], duplicated_genes)


#5b# Expression-based doublet detection with Scrublet, per sample
# Off by default, not dead code. Switch on when composition claims rest on exact
# cell counts. Left off here because 5c removes the doublets Scrublet cannot see
# and running both changes the denominator of every frequency. See docs/03_preprocessing.md.
RUN_SCRUBLET = False

if RUN_SCRUBLET:
    print("\nRunning Scrublet doublet detection (per sample)...")
    # scanpy moved scrublet out of .external around 1.10, so try both
    try:
        _scrublet = sc.pp.scrublet
    except AttributeError:
        _scrublet = sc.external.pp.scrublet

    for i, sample in enumerate(sampleNames):
        adata = adatas[i]
        n_before = int(adata.n_obs)

        # 10x multiplet rate is roughly 0.8% per 1,000 recovered cells, capped
        # so an unusually large lane cannot push the prior past a quarter
        edr = float(min(0.008 * n_before / 1000.0, 0.25))

        try:
            _scrublet(adata, expected_doublet_rate=edr, random_state=123)
            # predicted_doublet is NaN when Scrublet cannot auto-threshold
            pred = adata.obs["predicted_doublet"]
            adata.obs["predicted_doublet"] = pred.fillna(False).astype(bool)
            thr = adata.uns.get("scrublet", {}).get("threshold", np.nan)
        except Exception as e:
            # a failed sample must not abort the run, but it must be visible
            print(f"  Scrublet FAILED for {sample}: {e}\n  -> marking 0 doublets")
            adata.obs["predicted_doublet"] = False
            adata.obs["doublet_score"] = np.nan
            thr = np.nan

        n_doublets = int(adata.obs["predicted_doublet"].sum())
        dbl_rate = 100.0 * n_doublets / n_before if n_before else 0.0

        # retain singlets only
        adatas[i] = adata[~adata.obs["predicted_doublet"].values].copy()
        n_after = int(adatas[i].n_obs)

        print(f"  {sample}: {n_doublets}/{n_before} doublets "
              f"({dbl_rate:.2f}%), threshold={thr}, kept {n_after} singlets")

        qc_by_sample[sample].update({
            "scrublet_expected_rate": edr,
            "scrublet_threshold": float(thr) if thr == thr else np.nan,  # NaN-safe
            "n_doublets": n_doublets,
            "doublet_pct": dbl_rate,
            "n_cells_after_doublets": n_after,
        })

#----------------------------------------

#5c# TCR-based doublet detection (dual productive TRB), removing flagged barcodes
#
# Two distinct productive TRB chains on one barcode means two T cells in one droplet.
# This is a physical observation about the GEM, so it catches homotypic T-T doublets
# that Scrublet cannot. Alpha is not used (dual-alpha is common in real T cells) and a
# second beta must clear a UMI floor and a fraction of the dominant beta to rule out
# ambient contigs. Full reasoning and references: docs/03_preprocessing.md.

print("\nLocating vdj_t outputs...")
vdj_paths = []
for b in VDJ_BATCHES:
    vdj_paths += glob.glob(f"{VDJ_RAW_DIR}/{b}/cellranger_output/*/per_sample_outs/*/vdj_t/"f"filtered_contig_annotations.csv")

vdj_path_of = {}
for p in vdj_paths:
    hits = [s for t, s in VDJ_TOKEN_TO_SAMPLE.items() if t in p]
    if len(hits) == 1:
        vdj_path_of.setdefault(hits[0], []).append(p)
    else:
        print(f"  unresolved ({len(hits)} matches): {p}")

ambiguous = {s: v for s, v in vdj_path_of.items() if len(v) > 1}
if ambiguous:
    raise ValueError(f"more than one vdj_t file matched: {ambiguous}")

missing_vdj = sorted(set(sampleNames) - set(vdj_path_of))
print(f"  files found: {len(vdj_paths)} | samples resolved: {len(vdj_path_of)}/{len(sampleNames)}")
if missing_vdj:
    print(f"  WARNING no vdj_t, no TCR filtering applied: {missing_vdj}")


def flag_tcr_doublets(path):
    """Return (doublet_barcode_set, stats) for one filtered_contig_annotations.csv."""
    c = pd.read_csv(path)

    # Cell Ranger writes these as True/False, but the dtype varies by version
    def _flag(col):
        return c[col].astype(str).str.lower().isin(["true", "yes"])

    c = c[_flag("productive") & _flag("high_confidence") & _flag("is_cell") & _flag("full_length")]

    # collapse contigs that reassemble the same junction, so one chain is not counted twice
    c = c.drop_duplicates(subset=["barcode", "chain", "cdr3_nt"])

    beta = c[c["chain"] == "TRB"].sort_values(["barcode", "umis"], ascending=[True, False]).copy()
    beta["chain_rank"] = beta.groupby("barcode").cumcount()

    top    = beta.loc[beta["chain_rank"] == 0].set_index("barcode")["umis"]
    second = beta.loc[beta["chain_rank"] == 1].set_index("barcode")["umis"]

    supported = second[(second >= MIN_SECOND_BETA_UMIS)
                       & (second >= MIN_SECOND_BETA_FRAC * top.reindex(second.index))]

    n_tcr = c["barcode"].nunique()
    stats = {
        "n_tcr_barcodes":     int(n_tcr),
        "n_dual_beta_raw":    int(len(second)),
        "n_dual_beta_kept":   int(len(supported)),
        "n_dual_alpha":       int((c[c["chain"] == "TRA"].groupby("barcode").size() >= 2).sum()),
    }
    return set(supported.index), stats


print("\nApplying TCR doublet filter (per sample)...")
for i, sample in enumerate(sampleNames):
    adata = adatas[i]
    n_before = int(adata.n_obs)

    if sample not in vdj_path_of:
        qc_by_sample[sample].update({
            "n_tcr_barcodes": 0, "tcr_recovery_pct": np.nan, "n_dual_beta_raw": 0,
            "n_tcr_doublets": 0, "tcr_doublet_pct": 0.0, "n_cells_after_tcr": n_before,
        })
        adata.obs["has_tcr"] = False
        adata.obs["tcr_doublet"] = False
        print(f"  {sample}: no vdj_t, kept all {n_before} cells")
        continue

    dbl_bc, stats = flag_tcr_doublets(vdj_path_of[sample][0])

    # barcodes are still raw here; the sample suffix is appended in step 6
    contig_bc = pd.read_csv(vdj_path_of[sample][0], usecols=["barcode"])["barcode"].unique()
    adata.obs["has_tcr"] = adata.obs_names.isin(contig_bc)
    adata.obs["tcr_doublet"] = adata.obs_names.isin(dbl_bc)

    n_doublets = int(adata.obs["tcr_doublet"].sum())
    adatas[i] = adata[~adata.obs["tcr_doublet"].to_numpy()].copy()
    n_after = int(adatas[i].n_obs)

    tcr_pct = 100.0 * int(adata.obs["has_tcr"].sum()) / n_before if n_before else 0.0
    dbl_pct = 100.0 * n_doublets / n_before if n_before else 0.0

    print(f"  {sample}: TCR recovered on {tcr_pct:.1f}% of cells | "
          f"dual-beta {stats['n_dual_beta_raw']} raw -> {stats['n_dual_beta_kept']} supported | "
          f"dual-alpha {stats['n_dual_alpha']} (retained) | "
          f"dropped {n_doublets} ({dbl_pct:.2f}%), kept {n_after}")

    qc_by_sample[sample].update({
        "n_tcr_barcodes":     stats["n_tcr_barcodes"],
        "tcr_recovery_pct":   tcr_pct,
        "n_dual_beta_raw":    stats["n_dual_beta_raw"],
        "n_tcr_doublets":     n_doublets,
        "tcr_doublet_pct":    dbl_pct,
        "n_cells_after_tcr":  n_after,
    })

# sanity check: the observed rate should sit in the same range as the 10x multiplet expectation
# (~0.8% per 1,000 recovered cells), scaled down because only T-containing doublets are visible
_obs = sum(qc_by_sample[s].get("n_tcr_doublets", 0) for s in sampleNames)
_tot = sum(qc_by_sample[s]["n_cells_loaded"] for s in sampleNames)
_exp = sum(0.008 * qc_by_sample[s]["n_cells_loaded"] ** 2 / 1000.0 for s in sampleNames)
print(f"\nTCR doublets removed: {_obs:,} of {_tot:,} cells ({100.0 * _obs / _tot:.2f}%) | "
      f"10x-expected total multiplets across all types: ~{_exp:,.0f}")

#----------------------------------------

#6# Add an index to barcodes to indicate sample origin
for i, sample in enumerate(sampleNames):
    adatas[i].obs_names = adatas[i].obs_names + f"_{sample}"

#----------------------------------------

#7# Attach the design metadata from the sample sheet
# Joining rather than parsing: sample_origin is the key, every other design
# column is looked up. Nothing is inferred from the string form of the id.
meta_by_sample = samples_df.set_index("sample_origin")
meta_cols = ["patient_id", "timepoint", "timepoint_days",
             "category", "condition", "category_plot"]

for i, sample in enumerate(sampleNames):
    row = meta_by_sample.loc[sample]
    for col in meta_cols:
        adatas[i].obs[col] = row[col]

#8# Merge all AnnData objects into one
merged = ad.concat(adatas, merge="same")

# ordered categoricals so plots and model contrasts inherit the intended order
merged.obs["category"] = pd.Categorical(merged.obs["category"], categories=category_levels, ordered=True)
merged.obs["condition"] = pd.Categorical(merged.obs["condition"], categories=condition_levels, ordered=True)
merged.obs["category_plot"] = pd.Categorical(merged.obs["category_plot"], categories=category_plot_levels, ordered=True)

print(merged.obs["sample_origin"].value_counts())
print(merged.obs["patient_id"].value_counts())
print(merged.obs["timepoint"].value_counts())
print(merged.obs["category"].value_counts())
print(merged.obs["condition"].value_counts())
print(merged.obs["category_plot"].value_counts())

#--------------------------------------

#8b# Write the QC files (per-sample summary + UMI percentile table)
qc_df = pd.DataFrame([qc_by_sample[s] for s in sampleNames])
qc_df["retained_pct"] = 100.0 * qc_df["n_cells_after_tcr"] / qc_df["n_cells_loaded"]

# a total row so the aggregate is captured too
totals = {
    "sample_origin": "ALL",
    "n_cells_loaded": int(qc_df["n_cells_loaded"].sum()),
    "n_tcr_doublets": int(qc_df["n_tcr_doublets"].sum()),
    "n_cells_after_tcr": int(qc_df["n_cells_after_tcr"].sum()),
}
if RUN_SCRUBLET:
    totals["n_doublets"] = int(qc_df["n_doublets"].sum())
    totals["n_cells_after_doublets"] = int(qc_df["n_cells_after_doublets"].sum())
    totals["doublet_pct"] = 100.0 * totals["n_doublets"] / totals["n_cells_loaded"]
totals["tcr_doublet_pct"] = 100.0 * totals["n_tcr_doublets"] / totals["n_cells_loaded"]
totals["retained_pct"] = 100.0 * totals["n_cells_after_tcr"] / totals["n_cells_loaded"]
qc_df = pd.concat([qc_df, pd.DataFrame([totals])], ignore_index=True)

#QC workbook (02 appends to the same file)
qc_xlsx = f"{qc_csv_dir}/{PROJECT}_global_metacell_QC.xlsx"

def _write_sheets(path, sheets):
    """Create-or-append sheets to an .xlsx workbook (idempotent)."""
    if os.path.exists(path):
        with pd.ExcelWriter(path, engine="openpyxl", mode="a", if_sheet_exists="replace") as xw:
            for name, d in sheets.items():
                d.to_excel(xw, sheet_name=name[:31], index=False)
    else:
        with pd.ExcelWriter(path, engine="openpyxl", mode="w") as xw:
            for name, d in sheets.items():
                d.to_excel(xw, sheet_name=name[:31], index=False)

_write_sheets(qc_xlsx, {
    "preprocess_qc_summary": qc_df,
    "umi_percentiles": pd.DataFrame(umi_pct_rows),
})
print(f"\nWrote QC workbook - {qc_xlsx} (sheets: preprocess_qc_summary, umi_percentiles)")
print(qc_df.to_string(index=False))


#9# Save the final merged object as h5ad
out_path = f"{h5ad_dir}/{PROJECT}.merged_samples.h5ad"
merged.write_h5ad(out_path)
print(f"Saved merged object to {out_path}")
print("Successfully completed")
