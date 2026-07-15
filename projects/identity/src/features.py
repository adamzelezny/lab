from typing import Tuple, List
import numpy as np
import scanpy as sc
import anndata as ad
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


def split_train_test(
    adata: ad.AnnData,
    test_size: float = 0.2,
    random_state: int = 42
) -> Tuple[ad.AnnData, ad.AnnData]:
    indices = np.arange(adata.n_obs)
    train_idx, test_idx = train_test_split(
        indices,
        test_size=test_size,
        random_state=random_state,
        stratify=adata.obs["cell_type"]
    )
    adata_train = adata[train_idx].copy()
    adata_test = adata[test_idx].copy()
    return adata_train, adata_test


def select_highly_variable_genes(
    adata_train: ad.AnnData,
    n_top_genes: int = 3000
) -> List[str]:
    adata_tmp = adata_train.copy()
    sc.pp.highly_variable_genes(
        adata_tmp,
        n_top_genes=n_top_genes,
        flavor="seurat",
        subset=False
    )
    hvg_mask = adata_tmp.var["highly_variable"]
    hvg_list = adata_tmp.var_names[hvg_mask].tolist()
    return hvg_list


def extract_and_scale_features(
    adata_train: ad.AnnData,
    adata_test: ad.AnnData,
    selected_genes: List[str]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    X_train = adata_train[:, selected_genes].X
    X_test = adata_test[:, selected_genes].X
    
    if hasattr(X_train, "toarray"):
        X_train = X_train.toarray()
    if hasattr(X_test, "toarray"):
        X_test = X_test.toarray()
        
    y_train = adata_train.obs["cell_type"].values
    y_test = adata_test.obs["cell_type"].values
    
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    return X_train_scaled, y_train, X_test_scaled, y_test


def random_oversample(
    X: np.ndarray,
    y: np.ndarray,
    random_state: int = 42
) -> Tuple[np.ndarray, np.ndarray]:
    np.random.seed(random_state)
    unique_classes, class_counts = np.unique(y, return_counts=True)
    max_count = np.max(class_counts)
    
    X_resampled = []
    y_resampled = []
    
    for cls in unique_classes:
        indices = np.where(y == cls)[0]
        resampled_indices = np.random.choice(indices, size=max_count, replace=True)
        X_resampled.append(X[resampled_indices])
        y_resampled.append(y[resampled_indices])
        
    X_resampled = np.vstack(X_resampled)
    y_resampled = np.concatenate(y_resampled)
    
    shuffle_indices = np.random.permutation(len(y_resampled))
    X_resampled = X_resampled[shuffle_indices]
    y_resampled = y_resampled[shuffle_indices]
    
    return X_resampled, y_resampled
