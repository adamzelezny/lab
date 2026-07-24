#Author: Adam Zelezny
#Date: 2026-07-21

import sys
import os
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.stats
sys.path.append(os.path.join(os.path.dirname(__file__)))
from explore_gsea_class import compute_population_masks

def main():
    h5ad_path = os.path.join(os.path.dirname(__file__), '..', 'pik3ca_preprocessed.h5ad')
    adata = sc.read_h5ad(h5ad_path)
    tumor_mask, neuron_mask = compute_population_masks(adata)
    raw_adata = adata.raw.to_adata()
    shap_data = {'tumor_h1047r': [('Xist', 0.028019), ('EGFP', 0.022782), ('Slc1a2', 0.012173), ('Cdkn2a', 0.011225), ('Rpl12', 0.011087), ('Mt3', 0.01043), ('Egr1', 0.010393), ('Gm42418', 0.010266), ('Clu', 0.009523), ('1500015O10Rik', 0.009347), ('Cst3', 0.009039), ('Gap43', 0.008309), ('Serpina3n', 0.008087), ('Gpnmb', 0.00785), ('Lhfpl3', 0.006865)], 'tumor_c420r': [('Gm42418', 0.018638), ('Xist', 0.009505), ('Malat1', 0.007971), ('H2-Eb1', 0.007871), ('Mbp', 0.007735), ('EGFP', 0.00646), ('Mt3', 0.006192), ('Mt1', 0.005807), ('Fxyd1', 0.0055), ('Egr1', 0.00532), ('Tmsb4x', 0.005232), ('Slc1a2', 0.005183), ('Rpl12', 0.005139), ('Cst3', 0.005115), ('Cd74', 0.004937)], 'neuron_h1047r': [('EGFP', 0.021602), ('Gm42418', 0.013763), ('Met', 0.013184), ('Ttr', 0.012977), ('H2-Eb1', 0.011031), ('Egr1', 0.007917), ('Gpnmb', 0.006995), ('Meis2', 0.005534), ('Rpl12', 0.005484), ('H2-Aa', 0.005306), ('Xist', 0.005081), ('Cd74', 0.004995), ('Fos', 0.004758), ('Csrp1', 0.004365), ('H2-Ab1', 0.004084)], 'neuron_c420r': [('Gm42418', 0.015862), ('H2-Eb1', 0.010624), ('Cd74', 0.008232), ('Xist', 0.00722), ('Mbp', 0.006906), ('H2-Aa', 0.006294), ('Met', 0.005766), ('EGFP', 0.005568), ('Hist1h2bc', 0.00548), ('Malat1', 0.005404), ('Mog', 0.004898), ('Csrp1', 0.004739), ('Mal', 0.004698), ('Ttr', 0.004381), ('Thrsp', 0.003996)]}

    def calc_metrics(mask, items, g1='H1047R', g2='C420R'):
        m_g1 = mask & (adata.obs['variant'] == g1)
        m_g2 = mask & (adata.obs['variant'] == g2)
        sub1 = raw_adata[m_g1]
        sub2 = raw_adata[m_g2]
        lines = []
        for gene, shap_val in items:
            if gene not in raw_adata.var_names:
                continue
            g_idx = raw_adata.var_names.get_loc(gene)
            e1 = sub1.X[:, g_idx].toarray().flatten() if hasattr(sub1.X, 'toarray') else np.array(sub1.X[:, g_idx]).flatten()
            e2 = sub2.X[:, g_idx].toarray().flatten() if hasattr(sub2.X, 'toarray') else np.array(sub2.X[:, g_idx]).flatten()
            mean1 = float(np.mean(e1))
            mean2 = float(np.mean(e2))
            pct1 = float(np.mean(e1 > 0) * 100)
            pct2 = float(np.mean(e2 > 0) * 100)
            l2fc = float(np.log2((mean1 + 0.0001) / (mean2 + 0.0001)))
            pval = float(scipy.stats.ttest_ind(e1, e2, equal_var=False).pvalue)
            if np.isnan(pval):
                pval = 1.0
            lines.append(f'    {gene} & {mean1:.6f} & {mean2:.6f} & {pct1:.4f} & {pct2:.4f} & {l2fc:.6f} & {shap_val:.6f} \\\\')
        return '\n'.join(lines)
    print('--- TUMOR H1047R SHAP DRIVERS (H1047R vs C420R) ---')
    print(calc_metrics(tumor_mask, shap_data['tumor_h1047r'], g1='H1047R', g2='C420R'))
    print('\n--- TUMOR C420R SHAP DRIVERS (C420R vs H1047R) ---')
    print(calc_metrics(tumor_mask, shap_data['tumor_c420r'], g1='C420R', g2='H1047R'))
    print('\n--- NEURON H1047R SHAP DRIVERS (H1047R vs C420R) ---')
    print(calc_metrics(neuron_mask, shap_data['neuron_h1047r'], g1='H1047R', g2='C420R'))
    print('\n--- NEURON C420R SHAP DRIVERS (C420R vs H1047R) ---')
    print(calc_metrics(neuron_mask, shap_data['neuron_c420r'], g1='C420R', g2='H1047R'))
if __name__ == '__main__':
    main()
