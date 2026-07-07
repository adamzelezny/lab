import numpy as np
import pandas as pd

def create_mouse_brain_rna_dataset(
    n_cells=1000,
    random_state=42,
    save_path="synthetic_mouse_brain_cells.csv"
):
    np.random.seed(random_state)

    # Genes: some astrocyte markers, neuron markers, oligodendrocyte markers, microglia markers
    genes = [
        # Astrocyte markers
        "Gfap", "Aqp4", "Aldh1l1", "Slc1a3", "S100b",

        # Neuron markers
        "Snap25", "Rbfox3", "Map2", "Syt1", "Tubb3",

        # Oligodendrocyte markers
        "Mbp", "Mog", "Plp1", "Cnp", "Olig2",

        # Microglia markers
        "Cx3cr1", "P2ry12", "Tmem119", "Aif1", "Csf1r",

        # General / noisy genes
        "Actb", "Gapdh", "Rpl13a", "Malat1", "Jun"
    ]

    rows = []
    labels = []

    for i in range(n_cells):
        is_astrocyte = np.random.rand() < 0.35  # 35% astrocytes
        labels.append(int(is_astrocyte))

        # Determine sub-celltype if not astrocyte
        if not is_astrocyte:
            # 50% neurons, 35% oligodendrocytes, 15% microglia
            other_type = np.random.choice(["neuron", "oligodendrocyte", "microglia"], p=[0.5, 0.35, 0.15])
        else:
            other_type = "astrocyte"

        counts = {}
        # Cell-specific library size scaling factor (sequencing depth variation)
        depth_scale = np.random.lognormal(mean=0, sigma=0.25)

        for gene in genes:
            # Base expression level (lambda) depending on the actual cell type
            if gene in ["Gfap", "Aqp4", "Aldh1l1", "Slc1a3", "S100b"]:
                # Astrocyte markers
                lam = np.random.uniform(25, 70) if is_astrocyte else np.random.uniform(2, 12)

            elif gene in ["Snap25", "Rbfox3", "Map2", "Syt1", "Tubb3"]:
                # Neuron markers
                lam = np.random.uniform(25, 70) if other_type == "neuron" else np.random.uniform(2, 12)

            elif gene in ["Mbp", "Mog", "Plp1", "Cnp", "Olig2"]:
                # Oligodendrocyte markers
                lam = np.random.uniform(25, 70) if other_type == "oligodendrocyte" else np.random.uniform(2, 12)

            elif gene in ["Cx3cr1", "P2ry12", "Tmem119", "Aif1", "Csf1r"]:
                # Microglia markers
                lam = np.random.uniform(25, 70) if other_type == "microglia" else np.random.uniform(2, 12)

            else:
                # Housekeeping genes (expressed in all cells but with noise)
                lam = np.random.uniform(15, 45)

            # 1. Biological Overdispersion (Gamma-Poisson / Negative Binomial)
            # Use Gamma to draw a cell-gene specific rate, then Poisson
            dispersion = 1.5
            cell_gene_lam = np.random.gamma(shape=max(0.1, lam) / dispersion, scale=dispersion)
            
            # 2. Sequencing depth scale
            cell_gene_lam *= depth_scale

            count = np.random.poisson(cell_gene_lam)

            # 3. Dropout (zero-inflation) - probability of dropout decreases with expression level
            p_dropout = np.exp(-0.15 * count)
            if np.random.rand() < p_dropout:
                count = 0

            counts[gene] = count

        rows.append(counts)

    df = pd.DataFrame(rows)
    df["is_astrocyte"] = labels

    df.to_csv(save_path, index=False)
    return df


df = create_mouse_brain_rna_dataset()

print(df.head())
print(df["is_astrocyte"].value_counts())