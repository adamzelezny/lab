from pathlib import Path

import numpy as np
import scanpy as sc


# ============================================================
# CONFIG
# ============================================================

INPUT_FILE = Path(
    "projects/differences/outputs/normalized/combined_normalized.h5ad"
)

OUTPUT_FILE = Path(
    "projects/differences/outputs/log_transformed/combined_log1p.h5ad"
)

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)


# ============================================================
# LOAD
# ============================================================

print("\nLoading normalized dataset...")

adata = sc.read_h5ad(INPUT_FILE)

print(adata)


# ============================================================
# PRESERVE NORMALIZED MATRIX
# ============================================================

print("\nSaving normalized matrix into .layers['normalized']...")

adata.layers["normalized"] = adata.X.copy()

print("Done.")


# ============================================================
# BEFORE LOG TRANSFORM
# ============================================================

print("\n============================================================")
print("BEFORE LOG TRANSFORMATION")
print("============================================================")

values = adata.X

print(f"Minimum value : {values.min():.4f}")
print(f"Maximum value : {values.max():.4f}")

sample = values[:100, :100].toarray()

print("\nSample statistics:")
print(f"Mean : {sample.mean():.4f}")
print(f"Std  : {sample.std():.4f}")


# ============================================================
# LOG TRANSFORM
# ============================================================

print("\nApplying log1p transformation...")

sc.pp.log1p(adata)

print("Done.")


# ============================================================
# AFTER LOG TRANSFORM
# ============================================================

print("\n============================================================")
print("AFTER LOG TRANSFORMATION")
print("============================================================")

values = adata.X

print(f"Minimum value : {values.min():.4f}")
print(f"Maximum value : {values.max():.4f}")

sample = values[:100, :100].toarray()

print("\nSample statistics:")
print(f"Mean : {sample.mean():.4f}")
print(f"Std  : {sample.std():.4f}")


# ============================================================
# VERIFY LAYERS
# ============================================================

print("\n============================================================")
print("VERIFY DATA REPRESENTATIONS")
print("============================================================")

print("Layers present:")
print(list(adata.layers.keys()))

print("\nStored representations:")
print("- layers['counts']      -> raw integer counts")
print("- layers['normalized']  -> normalized counts")
print("- X                     -> log1p(normalized counts)")


# ============================================================
# SAVE
# ============================================================

print("\nSaving log-transformed dataset...")

adata.write_h5ad(OUTPUT_FILE)

print(f"\nSaved to:\n{OUTPUT_FILE}")