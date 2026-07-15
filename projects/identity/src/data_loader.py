import os
from typing import Tuple
import pandas as pd
import scanpy as sc
import anndata as ad


def load_raw_datasets(
    control_dir: str,
    metadata_path: str
) -> Tuple[ad.AnnData, pd.DataFrame]:
    adata = sc.read_10x_mtx(
        control_dir,
        var_names="gene_symbols",
        make_unique=True,
        cache=False
    )
    df_metadata = pd.read_csv(metadata_path)
    return adata, df_metadata


def align_metadata(
    adata: ad.AnnData,
    df_metadata: pd.DataFrame
) -> ad.AnnData:
    df_metadata = df_metadata.copy()
    df_metadata["barcode_stripped"] = df_metadata["barcode"].apply(
        lambda x: x.split("_")[0]
    )
    df_saline = df_metadata[
        df_metadata["barcode"].str.endswith("_8")
    ].drop_duplicates(subset=["barcode_stripped"])
    
    cell_type_map = dict(zip(df_saline["barcode_stripped"], df_saline["cellType"]))
    treatment_map = dict(zip(df_saline["barcode_stripped"], df_saline["treatment"]))
    
    adata.obs["cell_type"] = adata.obs_names.map(cell_type_map)
    adata.obs["treatment"] = adata.obs_names.map(treatment_map)
    adata_labeled = adata[adata.obs["cell_type"].notna()].copy()
    return adata_labeled


def filter_rare_cell_types(
    adata: ad.AnnData,
    min_cell_count: int = 10
) -> ad.AnnData:
    cell_type_counts = adata.obs["cell_type"].value_counts()
    keep_cell_types = cell_type_counts[cell_type_counts >= min_cell_count].index.tolist()
    adata_filtered = adata[adata.obs["cell_type"].isin(keep_cell_types)].copy()
    adata_filtered.obs["cell_type"] = adata_filtered.obs["cell_type"].astype("category")
    return adata_filtered
