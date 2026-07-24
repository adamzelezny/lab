#Author: Adam Zelezny
#Date: 2026-07-15

import scanpy as sc
import anndata as ad

def run_qc_filtering(adata: ad.AnnData, min_genes: int=500, max_genes: int=6000, max_mt_percent: float=15.0, min_counts: int=1000) -> ad.AnnData:
    adata_copy = adata.copy()
    adata_copy.var['mt'] = adata_copy.var_names.str.startswith('mt-')
    adata_copy.var['ribo'] = adata_copy.var_names.str.startswith(('Rps', 'Rpl'))
    sc.pp.calculate_qc_metrics(adata_copy, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
    keep_cells = (adata_copy.obs['n_genes_by_counts'] >= min_genes) & (adata_copy.obs['n_genes_by_counts'] <= max_genes) & (adata_copy.obs['pct_counts_mt'] <= max_mt_percent) & (adata_copy.obs['total_counts'] >= min_counts)
    adata_filtered = adata_copy[keep_cells].copy()
    return adata_filtered

def filter_rare_genes(adata: ad.AnnData, min_cells: int=3) -> ad.AnnData:
    adata_copy = adata.copy()
    sc.pp.filter_genes(adata_copy, min_cells=min_cells)
    return adata_copy

def normalize_and_log_transform(adata: ad.AnnData, target_sum: float=10000.0) -> ad.AnnData:
    adata_copy = adata.copy()
    adata_copy.layers['counts'] = adata_copy.X.copy()
    sc.pp.normalize_total(adata_copy, target_sum=target_sum)
    adata_copy.layers['normalized'] = adata_copy.X.copy()
    sc.pp.log1p(adata_copy)
    return adata_copy
