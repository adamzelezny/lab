from pathlib import Path

# pyrefly: ignore [missing-import]
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt


INPUT_FILE = Path("projects/differences/outputs/umap/combined_log1p_hvg_pca_umap_leiden.h5ad")

OUT_DIR = Path("projects/differences/outputs/markers")
PLOT_DIR = OUT_DIR / "plots"

OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)

CLUSTER_KEY = "leiden_0.5"
OUTPUT_FILE = OUT_DIR / "combined_with_cluster_markers.h5ad"


print("\nLoading UMAP + Leiden dataset...")
adata = sc.read_h5ad(INPUT_FILE)
print(adata)

print(f"\nUsing cluster key: {CLUSTER_KEY}")
print(adata.obs[CLUSTER_KEY].value_counts().sort_index())


print("\nFinding marker genes for each cluster...")
sc.tl.rank_genes_groups(
    adata,
    groupby=CLUSTER_KEY,
    method="wilcoxon",
    use_raw=False,
)

print("Marker gene analysis complete.")


print("\nSaving full marker table...")
markers = sc.get.rank_genes_groups_df(
    adata,
    group=None,
)

markers_file = OUT_DIR / f"{CLUSTER_KEY}_marker_genes_all.csv"
markers.to_csv(markers_file, index=False)
print(f"Saved: {markers_file}")


print("\nSaving top 20 markers per cluster...")
top20 = (
    markers
    .groupby("group")
    .head(20)
    .copy()
)

top20_file = OUT_DIR / f"{CLUSTER_KEY}_top20_marker_genes.csv"
top20.to_csv(top20_file, index=False)
print(f"Saved: {top20_file}")


print("\nTop 10 markers per cluster:")
for cluster in sorted(adata.obs[CLUSTER_KEY].unique(), key=lambda x: int(x)):
    cluster_markers = markers[markers["group"] == cluster].head(10)
    genes = cluster_markers["names"].tolist()
    print(f"Cluster {cluster}: {genes}")


print("\nCreating marker gene summary plot...")
sc.pl.rank_genes_groups(
    adata,
    n_genes=10,
    sharey=False,
    show=False,
)

plot_file = PLOT_DIR / f"{CLUSTER_KEY}_top_markers_plot.png"
plt.savefig(plot_file, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved: {plot_file}")


print("\nSaving updated AnnData object...")
adata.write_h5ad(OUTPUT_FILE)
print(f"Saved: {OUTPUT_FILE}")

print("\nDone.")
