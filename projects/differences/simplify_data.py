import os
import scanpy as sc
import pandas as pd
import numpy as np

CONTROL_PATH = r"C:\Users\adamz\Documents\lab\projects\differences\data\Control"
ACTIVITY_PATH = r"C:\Users\adamz\Documents\lab\projects\differences\data\Activity"
OUTPUT_CSV = r"C:\Users\adamz\Documents\lab\projects\differences\data\simplified_dataset.csv"

def main():
    print("Loading Control dataset...")
    control = sc.read_10x_mtx(CONTROL_PATH, var_names="gene_symbols", cache=True)
    
    print("Loading Activity dataset...")
    activity = sc.read_10x_mtx(ACTIVITY_PATH, var_names="gene_symbols", cache=True)
    
    # 1. Verify gene matching
    if not (control.var_names == activity.var_names).all():
        print("Warning: Gene lists do not match exactly. Aligning features...")
        # Find common genes
        common_genes = control.var_names.intersection(activity.var_names)
        control = control[:, common_genes]
        activity = activity[:, common_genes]
    
    genes = control.var_names.tolist()
    header = ["barcode", "condition"] + genes
    
    # WARNING: A dense matrix of 45,254 cells x 33,696 genes is huge (1.5 billion cells).
    # Converting the whole matrix to a dense pandas DataFrame all at once would require 
    # ~15+ GB of RAM and crash most systems.
    # 
    # SOLUTION: We stream and write the data to the CSV in chunks of 1000 cells at a time.
    print(f"Writing merged data to {OUTPUT_CSV} in chunks...")
    
    # Initialize the file with the header
    with open(OUTPUT_CSV, "w", encoding="utf-8") as f:
        f.write(",".join(header) + "\n")
        
    # Helper function to append cells in chunks
    def append_adata_to_csv(adata, condition, file_path, chunk_size=1000):
        n_cells = adata.n_obs
        for start in range(0, n_cells, chunk_size):
            end = min(start + chunk_size, n_cells)
            print(f"  Appending {condition} cells {start:,} to {end:,}...")
            
            # Slice the sparse matrix chunk and convert it to dense array
            chunk_dense = adata.X[start:end].toarray()
            chunk_barcodes = adata.obs_names[start:end]
            
            # Create a chunk DataFrame
            chunk_df = pd.DataFrame(chunk_dense, index=chunk_barcodes, columns=genes)
            chunk_df.insert(0, "condition", condition)
            chunk_df.index.name = "barcode"
            
            # Append to CSV (header=False since header is already written)
            chunk_df.to_csv(file_path, mode="a", header=False)
            
    append_adata_to_csv(control, "control", OUTPUT_CSV)
    append_adata_to_csv(activity, "activity", OUTPUT_CSV)
    
    print(f"\nSuccess! Merged dataset saved to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
