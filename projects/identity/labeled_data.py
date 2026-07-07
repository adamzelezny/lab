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

        counts = {}

        for gene in genes:
            if gene in ["Gfap", "Aqp4", "Aldh1l1", "Slc1a3", "S100b"]:
                # Astrocyte genes high in astrocytes, low otherwise
                lam = np.random.uniform(80, 250) if is_astrocyte else np.random.uniform(1, 20)

            elif gene in ["Snap25", "Rbfox3", "Map2", "Syt1", "Tubb3"]:
                # Neuron genes high in non-astrocytes
                lam = np.random.uniform(5, 25) if is_astrocyte else np.random.uniform(60, 200)

            elif gene in ["Mbp", "Mog", "Plp1", "Cnp", "Olig2"]:
                # Oligodendrocyte genes mostly non-astrocyte
                lam = np.random.uniform(2, 20) if is_astrocyte else np.random.uniform(40, 160)

            elif gene in ["Cx3cr1", "P2ry12", "Tmem119", "Aif1", "Csf1r"]:
                # Microglia genes mostly non-astrocyte
                lam = np.random.uniform(2, 15) if is_astrocyte else np.random.uniform(30, 140)

            else:
                # Housekeeping/noisy genes
                lam = np.random.uniform(20, 100)

            counts[gene] = np.random.poisson(lam)

        rows.append(counts)

    df = pd.DataFrame(rows)
    df["is_astrocyte"] = labels

    df.to_csv(save_path, index=False)
    return df


df = create_mouse_brain_rna_dataset()

print(df.head())
print(df["is_astrocyte"].value_counts())