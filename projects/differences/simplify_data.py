#Author: Adam Zelezny
#Date: 2026-07-15

import os
import scanpy as sc
import pandas as pd
import numpy as np
CONTROL_PATH = 'C:\\Users\\adamz\\Documents\\lab\\projects\\differences\\data\\Control'
ACTIVITY_PATH = 'C:\\Users\\adamz\\Documents\\lab\\projects\\differences\\data\\Activity'
OUTPUT_CSV = 'C:\\Users\\adamz\\Documents\\lab\\projects\\differences\\data\\simplified_dataset.csv'

def main():
    print('Loading Control dataset...')
    control = sc.read_10x_mtx(CONTROL_PATH, var_names='gene_symbols', cache=True)
    print('Loading Activity dataset...')
    activity = sc.read_10x_mtx(ACTIVITY_PATH, var_names='gene_symbols', cache=True)
    if not (control.var_names == activity.var_names).all():
        print('Warning: Gene lists do not match exactly. Aligning features...')
        common_genes = control.var_names.intersection(activity.var_names)
        control = control[:, common_genes]
        activity = activity[:, common_genes]
    genes = control.var_names.tolist()
    header = ['barcode', 'condition'] + genes
    print(f'Writing merged data to {OUTPUT_CSV} in chunks...')
    with open(OUTPUT_CSV, 'w', encoding='utf-8') as f:
        f.write(','.join(header) + '\n')

    def append_adata_to_csv(adata, condition, file_path, chunk_size=1000):
        n_cells = adata.n_obs
        for start in range(0, n_cells, chunk_size):
            end = min(start + chunk_size, n_cells)
            print(f'  Appending {condition} cells {start:,} to {end:,}...')
            chunk_dense = adata.X[start:end].toarray()
            chunk_barcodes = adata.obs_names[start:end]
            chunk_df = pd.DataFrame(chunk_dense, index=chunk_barcodes, columns=genes)
            chunk_df.insert(0, 'condition', condition)
            chunk_df.index.name = 'barcode'
            chunk_df.to_csv(file_path, mode='a', header=False)
    append_adata_to_csv(control, 'control', OUTPUT_CSV)
    append_adata_to_csv(activity, 'activity', OUTPUT_CSV)
    print(f'\nSuccess! Merged dataset saved to {OUTPUT_CSV}')
if __name__ == '__main__':
    main()
