import scanpy as sc
import anndata as ad


def run_qc_filtering(
    adata: ad.AnnData,
    min_genes: int = 500,
    max_genes: int = 6000,
    max_mt_percent: float = 15.0,
    min_counts: int = 1000
) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    
    keep_cells = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["pct_counts_mt"] <= max_mt_percent)
        & (adata.obs["total_counts"] >= min_counts)
    )
    
    adata_filtered = adata[keep_cells].copy()
    return adata_filtered


def filter_rare_genes(
    adata: ad.AnnData,
    min_cells: int = 3
) -> ad.AnnData:
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata


def normalize_and_log_transform(
    adata: ad.AnnData,
    target_sum: float = 10000.0
) -> ad.AnnData:
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    adata.layers["normalized"] = adata.X.copy()
    sc.pp.log1p(adata)
    return adata
