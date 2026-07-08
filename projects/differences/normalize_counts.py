from pathlib import Path

import numpy as np
import scanpy as sc


# ============================================================
# CONFIG
# ============================================================

INPUT_FILE = Path(
    "projects/differences/outputs/qc/combined_qc_filtered.h5ad"
)

OUTPUT_FILE = Path(
    "projects/differences/outputs/normalized/combined_normalized.h5ad"
)

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)


TARGET_SUM = 10000


# ============================================================
# LOAD
# ============================================================

print("\nLoading filtered dataset...")

adata = sc.read_h5ad(INPUT_FILE)

print(adata)


# ============================================================
# SAVE RAW COUNTS
# ============================================================

print("\nSaving original filtered counts into .layers['counts']...")

adata.layers["counts"] = adata.X.copy()

print("Done.")


# ============================================================
# BEFORE NORMALIZATION
# ============================================================

print("\n============================================================")
print("BEFORE NORMALIZATION")
print("============================================================")

counts_before = np.asarray(adata.X.sum(axis=1)).flatten()

print(f"Mean total counts : {counts_before.mean():.2f}")
print(f"Median            : {np.median(counts_before):.2f}")
print(f"Minimum           : {counts_before.min():.2f}")
print(f"Maximum           : {counts_before.max():.2f}")


# ============================================================
# NORMALIZE
# ============================================================

print("\nNormalizing every cell to total counts =", TARGET_SUM)

sc.pp.normalize_total(
    adata,
    target_sum=TARGET_SUM
)

print("Normalization complete.")


# ============================================================
# AFTER NORMALIZATION
# ============================================================

print("\n============================================================")
print("AFTER NORMALIZATION")
print("============================================================")

counts_after = np.asarray(adata.X.sum(axis=1)).flatten()

print(f"Mean total counts : {counts_after.mean():.2f}")
print(f"Median            : {np.median(counts_after):.2f}")
print(f"Minimum           : {counts_after.min():.2f}")
print(f"Maximum           : {counts_after.max():.2f}")

print("\nFirst 10 normalized totals:")

print(np.round(counts_after[:10], 2))


# ============================================================
# VERIFY RAW COUNTS LAYER
# ============================================================

print("\n============================================================")
print("VERIFY RAW COUNTS LAYER")
print("============================================================")

raw_counts = np.asarray(
    adata.layers["counts"].sum(axis=1)
).flatten()

print("Raw counts layer still contains:")

print(f"Mean total counts : {raw_counts.mean():.2f}")

print("\nCurrent X matrix contains:")

print(f"Mean total counts : {counts_after.mean():.2f}")


# ============================================================
# SAVE
# ============================================================

print("\nSaving normalized dataset...")

adata.write_h5ad(OUTPUT_FILE)

print(f"\nSaved to:\n{OUTPUT_FILE}")