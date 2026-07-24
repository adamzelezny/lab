#Author: Adam Zelezny
#Date: 2026-07-08

from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt
INPUT_FILE = Path('projects/differences/outputs/pca/combined_log1p_hvg_pca.h5ad')
OUT_DIR = Path('projects/differences/outputs/umap')
PLOT_DIR = OUT_DIR / 'plots'
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUT_DIR / 'combined_log1p_hvg_pca_umap_leiden.h5ad'
N_NEIGHBORS = 15
N_PCS = 50
LEIDEN_RESOLUTIONS = [0.2, 0.5, 1.0, 1.5]
print('\nLoading PCA dataset...')
adata = sc.read_h5ad(INPUT_FILE)
print(adata)
print('\nPCA shape:')
print(adata.obsm['X_pca'].shape)
print('\nBuilding nearest-neighbor graph...')
print(f'Using n_neighbors = {N_NEIGHBORS}')
print(f'Using n_pcs       = {N_PCS}')
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, use_rep='X_pca')
print('Neighbor graph complete.')
print('\nRunning Leiden clustering at multiple resolutions...')
for resolution in LEIDEN_RESOLUTIONS:
    key = f'leiden_{resolution}'
    print(f"Running Leiden resolution {resolution} -> adata.obs['{key}']")
    sc.tl.leiden(adata, resolution=resolution, key_added=key, flavor='igraph', directed=False, n_iterations=2)
    print(adata.obs[key].value_counts().sort_index())
print('\nComputing UMAP embedding...')
sc.tl.umap(adata)
print('UMAP complete.')
print("UMAP coordinates stored in adata.obsm['X_umap']")
print('UMAP shape:')
print(adata.obsm['X_umap'].shape)
print('\nSaving UMAP plots...')
plot_items = [('condition', 'umap_condition.png'), ('timepoint', 'umap_timepoint.png'), ('total_counts', 'umap_total_counts.png'), ('n_genes_by_counts', 'umap_genes_detected.png'), ('pct_counts_mt', 'umap_mito_percent.png'), ('pct_counts_ribo', 'umap_ribo_percent.png')]
for color_by, filename in plot_items:
    print(f'Plotting UMAP colored by {color_by}...')
    sc.pl.umap(adata, color=color_by, show=False, size=6, alpha=0.7)
    out_path = PLOT_DIR / filename
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')
for resolution in LEIDEN_RESOLUTIONS:
    key = f'leiden_{resolution}'
    filename = f'umap_{key}.png'
    print(f'Plotting UMAP colored by {key}...')
    sc.pl.umap(adata, color=key, legend_loc='on data', show=False, size=6, alpha=0.7)
    out_path = PLOT_DIR / filename
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')
print('\nSaving cluster composition tables...')
for resolution in LEIDEN_RESOLUTIONS:
    key = f'leiden_{resolution}'
    counts = adata.obs.groupby([key, 'condition']).size().unstack(fill_value=0)
    proportions = counts.div(counts.sum(axis=1), axis=0)
    counts_file = OUT_DIR / f'{key}_condition_counts.csv'
    proportions_file = OUT_DIR / f'{key}_condition_proportions.csv'
    counts.to_csv(counts_file)
    proportions.to_csv(proportions_file)
    print(f'Saved: {counts_file}')
    print(f'Saved: {proportions_file}')
print('\nSaving UMAP + Leiden dataset...')
adata.write_h5ad(OUTPUT_FILE)
print(f'\nSaved to:')
print(OUTPUT_FILE)
print('\nDone.')
