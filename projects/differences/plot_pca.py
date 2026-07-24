#Author: Adam Zelezny
#Date: 2026-07-08

from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt
INPUT_FILE = Path('projects/differences/outputs/pca/combined_log1p_hvg_pca.h5ad')
OUT_DIR = Path('projects/differences/outputs/pca/plots')
OUT_DIR.mkdir(parents=True, exist_ok=True)
print('\nLoading PCA dataset...')
adata = sc.read_h5ad(INPUT_FILE)
print(adata)
print('\nPCA coordinates stored in:')
print("adata.obsm['X_pca']")
print('\nPCA shape:')
print(adata.obsm['X_pca'].shape)
print('\nCreating PCA plots...')
plots = [('condition', 'pca_condition.png'), ('timepoint', 'pca_timepoint.png'), ('total_counts', 'pca_total_counts.png'), ('n_genes_by_counts', 'pca_genes_detected.png'), ('pct_counts_mt', 'pca_mito_percent.png'), ('pct_counts_ribo', 'pca_ribo_percent.png')]
for color_by, filename in plots:
    print(f'Plotting PCA colored by {color_by}...')
    sc.pl.pca(adata, color=color_by, components=['1,2'], show=False, size=8, alpha=0.6)
    out_path = OUT_DIR / filename
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')
print('\nCreating PCA variance ratio plot...')
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=False, show=False)
out_path = OUT_DIR / 'pca_variance_ratio.png'
plt.savefig(out_path, dpi=300, bbox_inches='tight')
plt.close()
print(f'Saved: {out_path}')
pc_pairs = [('1,2', 'pca_pc1_pc2_condition.png'), ('2,3', 'pca_pc2_pc3_condition.png'), ('3,4', 'pca_pc3_pc4_condition.png'), ('4,5', 'pca_pc4_pc5_condition.png')]
for components, filename in pc_pairs:
    print(f'Plotting PCA components {components} colored by condition...')
    sc.pl.pca(adata, color='condition', components=[components], show=False, size=8, alpha=0.6)
    out_path = OUT_DIR / filename
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')
print('\nDone.')
print(f'All PCA plots saved to: {OUT_DIR}')
