#Author: Adam Zelezny
#Date: 2026-07-15

import os
from typing import Tuple
import pandas as pd
import scanpy as sc
import anndata as ad

def load_raw_datasets(control_dir: str, metadata_path: str) -> Tuple[ad.AnnData, pd.DataFrame]:
    adata = sc.read_10x_mtx(control_dir, var_names='gene_symbols', make_unique=True, cache=False)
    df_metadata = pd.read_csv(metadata_path)
    return (adata, df_metadata)

def align_metadata(adata: ad.AnnData, df_metadata: pd.DataFrame, barcode_suffix: str='_8') -> ad.AnnData:
    df_metadata = df_metadata.copy()
    df_metadata['barcode_stripped'] = df_metadata['barcode'].apply(lambda x: x.split('_')[0])
    df_filtered = df_metadata[df_metadata['barcode'].str.endswith(barcode_suffix)].drop_duplicates(subset=['barcode_stripped'])
    cell_type_map = dict(zip(df_filtered['barcode_stripped'], df_filtered['cellType']))
    treatment_map = dict(zip(df_filtered['barcode_stripped'], df_filtered['treatment']))
    adata_copy = adata.copy()
    adata_copy.obs['cell_type'] = adata_copy.obs_names.map(cell_type_map)
    adata_copy.obs['treatment'] = adata_copy.obs_names.map(treatment_map)
    adata_labeled = adata_copy[adata_copy.obs['cell_type'].notna()].copy()
    return adata_labeled

def filter_rare_cell_types(adata: ad.AnnData, min_cell_count: int=10) -> ad.AnnData:
    cell_type_counts = adata.obs['cell_type'].value_counts()
    keep_cell_types = cell_type_counts[cell_type_counts >= min_cell_count].index.tolist()
    adata_filtered = adata[adata.obs['cell_type'].isin(keep_cell_types)].copy()
    adata_filtered.obs['cell_type'] = adata_filtered.obs['cell_type'].astype('category')
    return adata_filtered
