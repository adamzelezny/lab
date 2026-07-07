import scanpy as sc
import numpy as np

CONTROL_PATH = r"C:\Users\adamz\Documents\lab\projects\differences\data\Control"
ACTIVITY_PATH = r"C:\Users\adamz\Documents\lab\projects\differences\data\Activity"

print("Loading Control dataset...")
control = sc.read_10x_mtx(
    CONTROL_PATH,
    var_names="gene_symbols",
    cache=True
)

print("Loading Activity dataset...")
activity = sc.read_10x_mtx(
    ACTIVITY_PATH,
    var_names="gene_symbols",
    cache=True
)

print()


def summarize(name, adata):

    print("=" * 60)
    print(name)
    print("=" * 60)

    print(f"Cells : {adata.n_obs:,}")
    print(f"Genes : {adata.n_vars:,}")

    print()

    print("Expression matrix:")
    print(type(adata.X))

    print()

    nnz = adata.X.nnz
    total = adata.n_obs * adata.n_vars

    sparsity = 100 * (1 - nnz / total)

    print(f"Nonzero entries : {nnz:,}")
    print(f"Sparsity        : {sparsity:.2f}%")

    print()

    print("First 20 genes:")

    print(adata.var.index[:20].tolist())

    print()

    print("First 10 cell barcodes:")

    print(adata.obs.index[:10].tolist())

    print()

    print("Total counts per cell")

    counts = np.asarray(adata.X.sum(axis=1)).flatten()

    print(f"Mean   : {counts.mean():.1f}")
    print(f"Median : {np.median(counts):.1f}")
    print(f"Min    : {counts.min():.1f}")
    print(f"Max    : {counts.max():.1f}")

    print()

    print("Genes detected per cell")

    genes_detected = np.asarray((adata.X > 0).sum(axis=1)).flatten()

    print(f"Mean   : {genes_detected.mean():.1f}")
    print(f"Median : {np.median(genes_detected):.1f}")
    print(f"Min    : {genes_detected.min():.1f}")
    print(f"Max    : {genes_detected.max():.1f}")

    print()

    gene_counts = np.asarray(adata.X.sum(axis=0)).flatten()
    top = np.argsort(gene_counts)[::-1][:20]

    print("Top expressed genes:")
    for i in top:
        print(f"  {adata.var.index[i]}: {int(gene_counts[i]):,}")

    print()


summarize("CONTROL", control)

summarize("ACTIVITY", activity)