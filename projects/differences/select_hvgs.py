from pathlib import Path

import scanpy as sc


INPUT_FILE = Path("projects/differences/outputs/log_transformed/combined_log1p.h5ad")
OUTPUT_FILE = Path("projects/differences/outputs/hvg/combined_log1p_with_hvgs.h5ad")

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

N_TOP_GENES = 3000


print("\nLoading log-transformed dataset...")
adata = sc.read_h5ad(INPUT_FILE)
print(adata)

print("\nSelecting highly variable genes...")
print(f"Number of HVGs requested: {N_TOP_GENES:,}")

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=N_TOP_GENES,
    flavor="seurat",
    subset=False,
)

print("\nHVG selection complete.")

n_hvg = int(adata.var["highly_variable"].sum())
print(f"Highly variable genes selected: {n_hvg:,}")
print(f"Total genes retained in dataset: {adata.n_vars:,}")

print("\nTop 20 HVGs by normalized dispersion:")
top_hvgs = (
    adata.var[adata.var["highly_variable"]]
    .sort_values("dispersions_norm", ascending=False)
    .head(20)
)

print(top_hvgs[["means", "dispersions", "dispersions_norm"]])

print("\nImportant:")
print("- No genes were removed.")
print("- HVGs were marked in adata.var['highly_variable'].")
print("- Full gene matrix is still preserved.")

print(f"\nSaving dataset with HVG annotations to:\n{OUTPUT_FILE}")
adata.write_h5ad(OUTPUT_FILE)

print("\nDone.")