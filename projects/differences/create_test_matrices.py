import numpy as np
import pandas as pd

np.random.seed(42)

N_CELLS = 500
N_GENES = 500

#relevant gbm genes

important_genes = [
    "Sox2", "Gfap", "Olig2", "Pdgfra", "Mki67", "Top2a", "Pcna",
    "Vegfa", "Hif1a", "Egf", "Egfr", "Nes", "Aqp4", "Vim",
    "Fos", "Jun", "Egr1", "Arc", "Ntrk2", "Bdnf",
    "Cxcl10", "Ccl2", "Tnf", "Il6", "Cd68", "Aif1",
    "Snap25", "Syp", "Map2", "Rbfox3",
    "Mog", "Mbp", "Plp1", "Pecam1", "Kdr"
]

remaining_genes = [f"Gene_{i}" for i in range(N_GENES - len(important_genes))]
genes = important_genes + remaining_genes

#metadata

conditions = np.random.choice(["pre_activation", "activation"], size=N_CELLS)
mouse_ids = np.random.choice(["mouse_1", "mouse_2", "mouse_3", "mouse_4"], size=N_CELLS)

cell_types = np.random.choice(
    ["tumor_cell", "astrocyte_like", "oligodendrocyte_like", "microglia", "endothelial"],
    size=N_CELLS,
    p=[0.55, 0.18, 0.12, 0.10, 0.05]
)

metadata = pd.DataFrame({
    "cell_id": [f"cell_{i}" for i in range(N_CELLS)],
    "mouse_id": mouse_ids,
    "condition": conditions,
    "cell_type": cell_types
})

#rna counts simulation

counts = np.zeros((N_CELLS, N_GENES), dtype=int)

# gene specific average expression
base_means = np.random.gamma(shape=1.5, scale=2.0, size=N_GENES)

# highly expressed genes
base_means[:len(important_genes)] += np.random.uniform(1, 6, len(important_genes))

for g in range(N_GENES):
    mean_expression = base_means[g]

    # negative binomial like overdispersed counts
    dispersion = 1.5
    p = dispersion / (dispersion + mean_expression)

    counts[:, g] = np.random.negative_binomial(
        n=dispersion,
        p=p,
        size=N_CELLS
    )

#cell type specific expression

gene_idx = {gene: i for i, gene in enumerate(genes)}

def boost_gene(gene, mask, amount):
    idx = gene_idx[gene]
    counts[mask, idx] += np.random.poisson(amount, size=mask.sum())

tumor_mask = metadata["cell_type"] == "tumor_cell"
astro_mask = metadata["cell_type"] == "astrocyte_like"
oligo_mask = metadata["cell_type"] == "oligodendrocyte_like"
micro_mask = metadata["cell_type"] == "microglia"
endo_mask = metadata["cell_type"] == "endothelial"
activation_mask = metadata["condition"] == "activation"

for gene in ["Sox2", "Olig2", "Pdgfra", "Nes", "Vim"]:
    boost_gene(gene, tumor_mask, amount=6)

for gene in ["Mki67", "Top2a", "Pcna"]:
    boost_gene(gene, tumor_mask, amount=4)

for gene in ["Gfap", "Aqp4"]:
    boost_gene(gene, astro_mask, amount=10)

for gene in ["Mog", "Mbp", "Plp1"]:
    boost_gene(gene, oligo_mask, amount=10)

for gene in ["Cd68", "Aif1", "Ccl2", "Cxcl10", "Tnf"]:
    boost_gene(gene, micro_mask, amount=8)

for gene in ["Pecam1", "Kdr", "Vegfa"]:
    boost_gene(gene, endo_mask, amount=8)

for gene in ["Fos", "Jun", "Egr1", "Arc", "Bdnf"]:
    boost_gene(gene, activation_mask, amount=12)

for gene in ["Vegfa", "Hif1a", "Cxcl10", "Il6"]:
    boost_gene(gene, activation_mask, amount=5)

for gene in ["Mki67", "Top2a", "Pcna"]:
    idx = gene_idx[gene]
    reduction = np.random.poisson(2, size=activation_mask.sum())
    counts[activation_mask, idx] = np.maximum(
        0,
        counts[activation_mask, idx] - reduction
    )

dropout_prob = np.exp(-counts / 4)
dropout = np.random.random(counts.shape) < dropout_prob * 0.25
counts[dropout] = 0

count_df = pd.DataFrame(
    counts,
    index=metadata["cell_id"],
    columns=genes
)

count_df.to_csv("synthetic_glioblastoma_counts_500x500.csv")
metadata.to_csv("synthetic_glioblastoma_metadata.csv", index=False)

print("Created:")
print("synthetic_glioblastoma_counts_500x500.csv")
print("synthetic_glioblastoma_metadata.csv")
print()
print("Counts shape:", count_df.shape)
print("Metadata shape:", metadata.shape)
