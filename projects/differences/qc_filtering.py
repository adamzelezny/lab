from pathlib import Path

import scanpy as sc


# ============================================================
# CONFIG
# ============================================================

INPUT_FILE = Path("projects/differences/outputs/qc/combined_labeled_with_qc_metrics.h5ad")
OUT_DIR = Path("projects/differences/outputs/qc")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FILE = OUT_DIR / "combined_qc_filtered.h5ad"


# ============================================================
# FILTER THRESHOLDS
# ============================================================
# Based on your QC inspection:
# - n_genes < 500 removes mostly low-complexity cells
# - n_genes > 6000 removes possible doublets / extreme outliers
# - pct_counts_mt > 15 removes likely damaged/dying cells
# - total_counts < 1000 removes low-depth cells

MIN_GENES = 500
MAX_GENES = 6000
MAX_MT_PERCENT = 15
MIN_COUNTS = 1000


# ============================================================
# LOAD DATA
# ============================================================

print("\nLoading QC-annotated dataset...")
adata = sc.read_h5ad(INPUT_FILE)

print(adata)

starting_cells = adata.n_obs
starting_genes = adata.n_vars

print(f"\nStarting cells: {starting_cells:,}")
print(f"Starting genes: {starting_genes:,}")


# ============================================================
# DEFINE FILTERS
# ============================================================

print("\nApplying QC thresholds:")
print(f"Keep cells with n_genes_by_counts >= {MIN_GENES}")
print(f"Keep cells with n_genes_by_counts <= {MAX_GENES}")
print(f"Keep cells with pct_counts_mt <= {MAX_MT_PERCENT}")
print(f"Keep cells with total_counts >= {MIN_COUNTS}")

keep_cells = (
    (adata.obs["n_genes_by_counts"] >= MIN_GENES)
    & (adata.obs["n_genes_by_counts"] <= MAX_GENES)
    & (adata.obs["pct_counts_mt"] <= MAX_MT_PERCENT)
    & (adata.obs["total_counts"] >= MIN_COUNTS)
)


# ============================================================
# REPORT FILTERING EFFECT
# ============================================================

removed_cells = (~keep_cells).sum()
kept_cells = keep_cells.sum()

print("\nFiltering summary:")
print(f"Cells kept    : {kept_cells:,}")
print(f"Cells removed : {removed_cells:,}")
print(f"Percent removed: {100 * removed_cells / starting_cells:.2f}%")

print("\nCells removed by individual criterion:")
print(f"n_genes < {MIN_GENES}: {(adata.obs['n_genes_by_counts'] < MIN_GENES).sum():,}")
print(f"n_genes > {MAX_GENES}: {(adata.obs['n_genes_by_counts'] > MAX_GENES).sum():,}")
print(f"pct_mt > {MAX_MT_PERCENT}: {(adata.obs['pct_counts_mt'] > MAX_MT_PERCENT).sum():,}")
print(f"total_counts < {MIN_COUNTS}: {(adata.obs['total_counts'] < MIN_COUNTS).sum():,}")

print("\nCells kept by condition:")
print(adata.obs.loc[keep_cells, "condition"].value_counts())

print("\nCells removed by condition:")
print(adata.obs.loc[~keep_cells, "condition"].value_counts())


# ============================================================
# APPLY FILTER
# ============================================================

adata_filtered = adata[keep_cells].copy()

print("\nFiltered dataset:")
print(adata_filtered)


# ============================================================
# OPTIONAL: FILTER GENES EXPRESSED IN VERY FEW CELLS
# ============================================================
# This removes genes detected in fewer than 3 cells after cell QC.

print("\nFiltering genes expressed in fewer than 3 cells...")

before_gene_filter = adata_filtered.n_vars
sc.pp.filter_genes(adata_filtered, min_cells=3)
after_gene_filter = adata_filtered.n_vars

print(f"Genes before: {before_gene_filter:,}")
print(f"Genes after : {after_gene_filter:,}")
print(f"Genes removed: {before_gene_filter - after_gene_filter:,}")


# ============================================================
# SAVE FILTERED DATA
# ============================================================

print(f"\nSaving filtered dataset to: {OUTPUT_FILE}")
adata_filtered.write_h5ad(OUTPUT_FILE)

print("\nDone.")
print("This file is now ready for normalization and log transformation.")