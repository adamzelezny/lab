#Author: Adam Zelezny
#Date: 2026-07-07

import numpy as np
import pandas as pd

def create_mouse_brain_rna_dataset(n_cells=1000, random_state=42, save_path='synthetic_mouse_brain_cells.csv'):
    np.random.seed(random_state)
    genes = ['Gfap', 'Aqp4', 'Aldh1l1', 'Slc1a3', 'S100b', 'Snap25', 'Rbfox3', 'Map2', 'Syt1', 'Tubb3', 'Mbp', 'Mog', 'Plp1', 'Cnp', 'Olig2', 'Cx3cr1', 'P2ry12', 'Tmem119', 'Aif1', 'Csf1r', 'Actb', 'Gapdh', 'Rpl13a', 'Malat1', 'Jun']
    rows = []
    labels = []
    for i in range(n_cells):
        is_astrocyte = np.random.rand() < 0.35
        labels.append(int(is_astrocyte))
        if not is_astrocyte:
            other_type = np.random.choice(['neuron', 'oligodendrocyte', 'microglia'], p=[0.5, 0.35, 0.15])
        else:
            other_type = 'astrocyte'
        counts = {}
        depth_scale = np.random.lognormal(mean=0, sigma=0.25)
        for gene in genes:
            if gene in ['Gfap', 'Aqp4', 'Aldh1l1', 'Slc1a3', 'S100b']:
                lam = np.random.uniform(25, 70) if is_astrocyte else np.random.uniform(2, 12)
            elif gene in ['Snap25', 'Rbfox3', 'Map2', 'Syt1', 'Tubb3']:
                lam = np.random.uniform(25, 70) if other_type == 'neuron' else np.random.uniform(2, 12)
            elif gene in ['Mbp', 'Mog', 'Plp1', 'Cnp', 'Olig2']:
                lam = np.random.uniform(25, 70) if other_type == 'oligodendrocyte' else np.random.uniform(2, 12)
            elif gene in ['Cx3cr1', 'P2ry12', 'Tmem119', 'Aif1', 'Csf1r']:
                lam = np.random.uniform(25, 70) if other_type == 'microglia' else np.random.uniform(2, 12)
            else:
                lam = np.random.uniform(15, 45)
            dispersion = 1.5
            cell_gene_lam = np.random.gamma(shape=max(0.1, lam) / dispersion, scale=dispersion)
            cell_gene_lam *= depth_scale
            count = np.random.poisson(cell_gene_lam)
            p_dropout = np.exp(-0.15 * count)
            if np.random.rand() < p_dropout:
                count = 0
            counts[gene] = count
        rows.append(counts)
    df = pd.DataFrame(rows)
    df['is_astrocyte'] = labels
    df.to_csv(save_path, index=False)
    return df
df = create_mouse_brain_rna_dataset()
print(df.head())
print(df['is_astrocyte'].value_counts())
