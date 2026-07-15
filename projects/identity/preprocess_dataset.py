import os
import sys
import numpy as np

project_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(project_dir)

from src.data_loader import load_raw_datasets, align_metadata, filter_rare_cell_types
from src.preprocessing import run_qc_filtering, filter_rare_genes, normalize_and_log_transform
from src.features import split_train_test, select_highly_variable_genes, extract_and_scale_features


def main():
    control_dir = os.path.join(project_dir, "data", "control")
    metadata_path = os.path.join(project_dir, "data", "p50_cell_metadata.csv")
    
    output_dir = os.path.join(project_dir, "data", "processed")
    os.makedirs(output_dir, exist_ok=True)
    
    adata, df_meta = load_raw_datasets(control_dir, metadata_path)
    adata_labeled = align_metadata(adata, df_meta)
    adata_filtered_types = filter_rare_cell_types(adata_labeled, min_cell_count=10)
    
    adata_qc = run_qc_filtering(
        adata_filtered_types,
        min_genes=500,
        max_genes=6000,
        max_mt_percent=15.0,
        min_counts=1000
    )
    adata_qc = filter_rare_genes(adata_qc, min_cells=3)
    adata_preprocessed = normalize_and_log_transform(adata_qc, target_sum=10000.0)
    
    adata_train, adata_test = split_train_test(adata_preprocessed, test_size=0.2, random_state=42)
    hvg_list = select_highly_variable_genes(adata_train, n_top_genes=3000)
    
    X_train, y_train, X_test, y_test = extract_and_scale_features(
        adata_train,
        adata_test,
        selected_genes=hvg_list
    )
    
    train_out_path = os.path.join(output_dir, "train_features.npz")
    test_out_path = os.path.join(output_dir, "test_features.npz")
    hvg_out_path = os.path.join(output_dir, "selected_hvg_genes.txt")
    
    np.savez_compressed(
        train_out_path,
        X=X_train,
        y=y_train
    )
    np.savez_compressed(
        test_out_path,
        X=X_test,
        y=y_test
    )
    
    with open(hvg_out_path, "w") as f:
        f.write("\n".join(hvg_list))


if __name__ == "__main__":
    main()
