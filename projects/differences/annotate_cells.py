#Author: Adam Zelezny
#Date: 2026-07-09

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
INPUT_FILE = Path('projects/differences/outputs/markers/combined_with_cluster_markers.h5ad')
OUT_DIR = Path('projects/differences/outputs/annotated')
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUT_DIR / 'combined_annotated.h5ad'
PROPORTIONS_PLOT_FILE = OUT_DIR / 'cell_type_proportions.png'
PROPORTIONS_CSV_FILE = OUT_DIR / 'cell_type_proportions.csv'
cluster_to_cell_type = {'0': 'T cells', '1': 'T cells (lymphoid)', '2': 'MHC-II+ Microglia/Macrophages', '3': 'Border-associated Macrophages', '4': 'MHC-II+ Microglia/Macrophages', '5': 'Endothelial cells', '6': 'Ependymal/Choroid plexus', '7': 'B cells', '8': 'Microglia/Macrophages', '9': 'Proliferating cells', '10': 'Dendritic cells', '11': 'Mature oligodendrocytes', '12': 'OPC-like cells', '13': 'Myelinating oligodendrocytes', '14': 'Pericytes/SMCs', '15': 'Fibroblasts', '16': 'Activated microglia', '17': 'Ependymal/Choroid plexus', '18': 'Homeostatic microglia', '19': 'Immature neurons', '20': 'Mature oligodendrocytes', '21': 'Astrocyte-like cells', '22': 'Ependymal cells', '23': 'Astrocyte/Microglia mixed', '24': 'Microglia (ribosomal-high)'}
print('Loading dataset...')
adata = sc.read_h5ad(INPUT_FILE)
print(adata)
print('\nAdding cell type annotations...')
adata.obs['cell_type'] = adata.obs['leiden_0.5'].map(cluster_to_cell_type)
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
print('\nCell type counts:')
print(adata.obs['cell_type'].value_counts())
print('\nCalculating condition proportions...')
counts_table = adata.obs.groupby(['cell_type', 'condition'], observed=False).size().unstack(fill_value=0)
proportions_table = counts_table.div(counts_table.sum(axis=0), axis=1)
print('\nCell type proportions by condition:')
print(proportions_table.round(4))
proportions_table.to_csv(PROPORTIONS_CSV_FILE)
print(f'Saved proportions table to: {PROPORTIONS_CSV_FILE}')
print('\nPlotting cell type proportions...')
fig, ax = plt.subplots(figsize=(10, 6))
proportions_table.T.plot(kind='bar', stacked=True, ax=ax, colormap='tab20', width=0.6)
ax.set_title('Cell Type Composition by Condition', fontsize=14, pad=15)
ax.set_ylabel('Proportion of total cells', fontsize=12)
ax.set_xlabel('Condition', fontsize=12)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Cell Types')
plt.tight_layout()
plt.savefig(PROPORTIONS_PLOT_FILE, dpi=300, bbox_inches='tight')
plt.close()
print(f'Saved proportions plot to: {PROPORTIONS_PLOT_FILE}')
print('\nSaving annotated AnnData dataset...')
adata.write_h5ad(OUTPUT_FILE)
print(f'Saved: {OUTPUT_FILE}')
print('\nDone.')
