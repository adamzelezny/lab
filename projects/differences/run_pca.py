#Author: Adam Zelezny
#Date: 2026-07-08

from pathlib import Path
import numpy as np
import scanpy as sc
import pandas as pd
INPUT_FILE = Path('projects/differences/outputs/hvg/combined_log1p_with_hvgs.h5ad')
OUT_DIR = Path('projects/differences/outputs/pca')
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUT_DIR / 'combined_log1p_hvg_pca.h5ad'
PCA_VARIANCE_FILE = OUT_DIR / 'pca_variance_ratio.csv'
N_PCS = 50
print('\nLoading HVG-annotated dataset...')
adata = sc.read_h5ad(INPUT_FILE)
print(adata)
n_hvgs = int(adata.var['highly_variable'].sum())
print(f'\nHighly variable genes available: {n_hvgs:,}')
print('\nCreating HVG-only view for PCA...')
adata_hvg = adata[:, adata.var['highly_variable']].copy()
print(adata_hvg)
print('\nScaling HVG expression values...')
sc.pp.scale(adata_hvg, max_value=10)
print('Scaling complete.')
print(f'\nRunning PCA with {N_PCS} components...')
sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver='arpack')
print('PCA complete.')
print('\nCopying PCA results back into full AnnData object...')
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
pcs_loadings = np.zeros((adata.n_vars, N_PCS))
highly_variable_mask = adata.var['highly_variable']
pcs_loadings[highly_variable_mask, :] = adata_hvg.varm['PCs']
adata.varm['PCs'] = pcs_loadings
adata.uns['pca'] = adata_hvg.uns['pca']
print("PCA stored in adata.obsm['X_pca'].")
print('\nPCA variance explained:')
variance_ratio = adata.uns['pca']['variance_ratio']
pca_df = pd.DataFrame({'PC': [f'PC{i + 1}' for i in range(len(variance_ratio))], 'variance_ratio': variance_ratio, 'cumulative_variance': variance_ratio.cumsum()})
print(pca_df.head(20).to_string(index=False))
pca_df.to_csv(PCA_VARIANCE_FILE, index=False)
print(f'\nSaved PCA variance table to:')
print(PCA_VARIANCE_FILE)
print('\nSaving PCA dataset...')
adata.write_h5ad(OUTPUT_FILE)
print(f'\nSaved to:')
print(OUTPUT_FILE)
print('\nDone.')
