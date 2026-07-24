#Author: Adam Zelezny
#Date: 2026-07-14

import os
import glob
import scanpy as sc
import numpy as np

def find_samples(data_dir: str) -> list[str]:
    pattern = os.path.join(data_dir, '*_matrix.mtx.gz')
    files = glob.glob(pattern)
    samples = []
    for f in files:
        basename = os.path.basename(f)
        sample = basename.replace('_matrix.mtx.gz', '')
        samples.append(sample)
    return sorted(samples)

def load_and_filter_sample(data_dir: str, sample_name: str, min_counts: int=500) -> sc.AnnData:
    adata = sc.read_10x_mtx(data_dir, prefix=f'{sample_name}_', cache=False)
    counts = np.array(adata.X.sum(axis=1)).flatten()
    adata = adata[counts >= min_counts].copy()
    adata.obs_names = f'{sample_name}_' + adata.obs_names
    adata.obs['sample'] = sample_name
    parts = sample_name.split('-')
    adata.obs['variant'] = parts[0]
    adata.obs['replicate'] = parts[1]
    print(f'Loaded {sample_name}: {adata.n_obs} cells')
    return adata

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    data_dir = os.path.join(project_dir, 'data')
    output_path = os.path.join(project_dir, 'pik3ca_merged.h5ad')
    samples = find_samples(data_dir)
    adatas = []
    for sample in samples:
        adata = load_and_filter_sample(data_dir, sample, min_counts=500)
        adatas.append(adata)
    adata_merged = sc.concat(adatas, join='inner')
    first_adata = adatas[0]
    adata_merged.var = first_adata.var.copy()
    print(f'Merged: {adata_merged.n_obs} cells x {adata_merged.n_vars} genes')
    print(f'Saving merged dataset to {output_path}')
    adata_merged.write_h5ad(output_path, compression='gzip')
if __name__ == '__main__':
    main()
