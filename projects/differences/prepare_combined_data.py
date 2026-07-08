from pathlib import Path
import scanpy as sc
import anndata as ad
import numpy as np


# ============================================================
# CONFIGURATION
# ============================================================

DATA_DIRECTORY = Path("projects/differences/data")
CONTROL_DATA_PATH = DATA_DIRECTORY / "Control"
ACTIVITY_DATA_PATH = DATA_DIRECTORY / "Activity"
COMBINED_OUTPUT_PATH = DATA_DIRECTORY / "combined_labeled_data.h5ad"


# ============================================================
# FUNCTIONS
# ============================================================

def load_dataset_from_10x_directory(directory_path: Path) -> ad.AnnData:
    """Loads a 10x Genomics MTX dataset into an AnnData object."""
    print(f"Loading 10x dataset from: {directory_path}")
    return sc.read_10x_mtx(
        directory_path,
        var_names="gene_symbols",
        cache=False
    )


def annotate_cell_metadata(
    dataset: ad.AnnData,
    condition: str,
    mouse_id: str,
    timepoint: str,
    source_folder: str
) -> None:
    """Adds biological and experimental metadata columns to the dataset."""
    print(f"Adding metadata for condition '{condition}' from '{source_folder}'...")
    dataset.obs["condition"] = condition
    dataset.obs["mouse_id"] = mouse_id
    dataset.obs["timepoint"] = timepoint
    dataset.obs["source_folder"] = source_folder


def make_barcodes_unique(dataset: ad.AnnData, prefix: str) -> None:
    """Prepends a prefix to the cell barcodes (obs_names) to make them unique."""
    print(f"Prepending prefix '{prefix}' to cell barcodes...")
    dataset.obs_names = [f"{prefix}_{barcode}" for barcode in dataset.obs_names]


def check_gene_compatibility(control_dataset: ad.AnnData, activity_dataset: ad.AnnData) -> None:
    """Verifies that both datasets have matching gene counts and gene ordering."""
    print("Verifying gene compatibility between control and activity datasets...")
    
    same_gene_count = control_dataset.n_vars == activity_dataset.n_vars
    same_gene_order = np.array_equal(control_dataset.var_names, activity_dataset.var_names)

    print(f"Same gene count: {same_gene_count}")
    print(f"Same gene order: {same_gene_order}")

    if not same_gene_order:
        raise ValueError(
            "Datasets are incompatible: Control and Activity datasets must have the exact same genes in the same order."
        )
    print("Gene compatibility check passed.")


def print_dataset_summary(dataset_name: str, dataset: ad.AnnData) -> None:
    """Prints a clear and simple summary of the dataset size, sparsity, and metadata."""
    print(f"\n--- Summary for dataset: {dataset_name} ---")
    
    cell_count = dataset.n_obs
    gene_count = dataset.n_vars
    print(f"Number of cells: {cell_count:,}")
    print(f"Number of genes: {gene_count:,}")
    print(f"Expression matrix type: {type(dataset.X)}")

    nonzero_entries_count = dataset.X.nnz
    total_matrix_elements = cell_count * gene_count
    sparsity_percentage = 100 * (1 - nonzero_entries_count / total_matrix_elements)

    print(f"Non-zero expression count: {nonzero_entries_count:,}")
    print(f"Matrix sparsity: {sparsity_percentage:.2f}%")

    print("First 5 cell barcodes:", list(dataset.obs_names[:5]))
    print("First 5 gene symbols:", list(dataset.var_names[:5]))
    print("Cell metadata columns:", list(dataset.obs.columns))
    print("Gene metadata columns:", list(dataset.var.columns))


def concatenate_datasets(control_dataset: ad.AnnData, activity_dataset: ad.AnnData) -> ad.AnnData:
    """Concatenates the control and activity AnnData objects along the cells axis."""
    print("Combining control and activity datasets...")
    combined_dataset = ad.concat(
        [control_dataset, activity_dataset],
        join="inner",
        label="dataset",
        keys=["control", "activity"],
        index_unique=None
    )
    combined_dataset.var = control_dataset.var.copy()
    return combined_dataset


def check_raw_counts_format(dataset: ad.AnnData) -> None:
    """Checks a small slice of the expression matrix to verify if they are raw integer counts."""
    print("\n--- Verifying raw count format (expecting integers) ---")
    
    subset_expression_matrix = dataset.X[:100, :100].toarray()
    minimum_count_value = subset_expression_matrix.min()
    maximum_count_value = subset_expression_matrix.max()
    unique_subset_values = np.unique(subset_expression_matrix[:10, :10])

    print(f"Minimum count value in slice: {minimum_count_value}")
    print(f"Maximum count value in slice: {maximum_count_value}")
    print(f"Unique values in first 10x10 submatrix: {unique_subset_values}")

    if np.all(subset_expression_matrix == subset_expression_matrix.astype(int)):
        print("Interpretation: Matrix contains raw integer counts (correct format).")
    else:
        print("Interpretation: Matrix contains decimal numbers (normalization or scaling already applied).")


def save_anndata_to_file(dataset: ad.AnnData, file_path: Path) -> None:
    """Saves the AnnData object to the specified file path in h5ad format."""
    print(f"Saving combined dataset to: {file_path}")
    dataset.write_h5ad(file_path)
    print("File saved successfully.")


def main() -> None:
    # 1. Load Data
    control_dataset = load_dataset_from_10x_directory(CONTROL_DATA_PATH)
    activity_dataset = load_dataset_from_10x_directory(ACTIVITY_DATA_PATH)

    # 2. Add Metadata
    annotate_cell_metadata(
        dataset=control_dataset,
        condition="control",
        mouse_id="mouse_A",
        timepoint="pre_activation",
        source_folder="Control"
    )
    annotate_cell_metadata(
        dataset=activity_dataset,
        condition="activity",
        mouse_id="mouse_A",
        timepoint="post_activation",
        source_folder="Activity"
    )

    make_barcodes_unique(control_dataset, prefix="control")
    make_barcodes_unique(activity_dataset, prefix="activity")

    print_dataset_summary("Control Dataset", control_dataset)
    print_dataset_summary("Activity Dataset", activity_dataset)

    check_gene_compatibility(control_dataset, activity_dataset)

    combined_dataset = concatenate_datasets(control_dataset, activity_dataset)

    print_dataset_summary("Combined Dataset", combined_dataset)
    
    print("\n--- Combined Metadata Value Counts ---")
    print("Condition counts:")
    print(combined_dataset.obs["condition"].value_counts())
    print("\nMouse ID counts:")
    print(combined_dataset.obs["mouse_id"].value_counts())
    print("\nTimepoint counts:")
    print(combined_dataset.obs["timepoint"].value_counts())
    print("\nSource folder counts:")
    print(combined_dataset.obs["source_folder"].value_counts())

    print("\nFirst 10 rows of combined cell metadata:")
    print(combined_dataset.obs.head(10))

    print("\nFirst 5 rows of combined gene metadata:")
    print(combined_dataset.var.head(5))

    # 8. Check Count Matrix Format
    check_raw_counts_format(combined_dataset)

    # 9. Save Combined Dataset
    save_anndata_to_file(combined_dataset, COMBINED_OUTPUT_PATH)


if __name__ == "__main__":
    main()