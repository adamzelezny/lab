#Author: Adam Zelezny
#Date: 2026-07-17

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.stats

def load_dataset(h5ad_path):
    return sc.read_h5ad(h5ad_path)

def compute_population_masks(adata):
    raw_adata = adata.raw.to_adata()
    gfap_idx = raw_adata.var_names.get_loc('Gfap')
    map2_idx = raw_adata.var_names.get_loc('Map2')
    gpc3_idx = raw_adata.var_names.get_loc('Gpc3')
    if hasattr(raw_adata.X, 'toarray'):
        gfap_expr = raw_adata.X[:, gfap_idx].toarray().flatten()
        map2_expr = raw_adata.X[:, map2_idx].toarray().flatten()
        gpc3_expr = raw_adata.X[:, gpc3_idx].toarray().flatten()
    else:
        gfap_expr = np.array(raw_adata.X[:, gfap_idx]).flatten()
        map2_expr = np.array(raw_adata.X[:, map2_idx]).flatten()
        gpc3_expr = np.array(raw_adata.X[:, gpc3_idx]).flatten()
    tumor_mask = (adata.obs['kmeans'] == '2') & ((gfap_expr > 0) | (gpc3_expr > 0)) & (map2_expr == 0)
    neuron_mask = (map2_expr > 0) & (gfap_expr == 0)
    return (tumor_mask, neuron_mask)

def select_comparison():
    print('\nSelect comparison:')
    print('1. H1047R vs R88Q (Control)')
    print('2. C420R vs R88Q (Control)')
    print('3. H1047R vs C420R (Direct)')
    choice = input('Enter choice (1-3): ').strip()
    if choice == '1':
        return ('H1047R', 'R88Q', '_H1047R_vs_R88Q_control')
    elif choice == '2':
        return ('C420R', 'R88Q', '_C420R_vs_R88Q_control')
    elif choice == '3':
        return ('H1047R', 'C420R', '_H1047R_vs_C420R_direct')
    else:
        raise ValueError('Invalid comparison selection')

def select_population():
    print('\nSelect cell population:')
    print('1. Tumor')
    print('2. Neurons')
    choice = input('Enter choice (1-2): ').strip()
    if choice == '1':
        return ('tumor', 'tumor')
    elif choice == '2':
        return ('neuron', 'neuron')
    else:
        raise ValueError('Invalid population selection')

def find_matching_pathway(df, query):
    df_matched = df[df['Term'].str.contains(query, case=False, na=False)]
    if df_matched.empty:
        return None
    return df_matched.iloc[0]

def compute_differential_metrics(adata, mask, group_a, group_b, genes):
    raw_adata = adata.raw.to_adata()[mask]
    results = []
    for gene in genes:
        if gene not in raw_adata.var_names:
            continue
        gene_idx = raw_adata.var_names.get_loc(gene)
        cells_a = (raw_adata.obs['variant'] == group_a).values
        cells_b = (raw_adata.obs['variant'] == group_b).values
        if hasattr(raw_adata.X, 'toarray'):
            expr_a = raw_adata.X[cells_a, gene_idx].toarray().flatten()
            expr_b = raw_adata.X[cells_b, gene_idx].toarray().flatten()
        else:
            expr_a = np.array(raw_adata.X[cells_a, gene_idx]).flatten()
            expr_b = np.array(raw_adata.X[cells_b, gene_idx]).flatten()
        mean_a = float(np.mean(expr_a))
        mean_b = float(np.mean(expr_b))
        pct_a = float(np.mean(expr_a > 0) * 100)
        pct_b = float(np.mean(expr_b > 0) * 100)
        log2_fc = float(np.log2(mean_a + 0.0001) - np.log2(mean_b + 0.0001))
        if len(expr_a) > 0 and len(expr_b) > 0:
            stat, p_val = scipy.stats.mannwhitneyu(expr_a, expr_b, alternative='two-sided')
        else:
            p_val = 1.0
        results.append({'Gene': gene, f'Mean_{group_a}': mean_a, f'Mean_{group_b}': mean_b, f'Pct_{group_a}': pct_a, f'Pct_{group_b}': pct_b, 'Log2FoldChange': log2_fc, 'PValue': p_val})
    return pd.DataFrame(results)

def show_top_pathways(df_gsea, count=10):
    df_gsea = df_gsea.copy()
    df_gsea['Abs_NES'] = df_gsea['NES'].abs()
    df_top = df_gsea.sort_values(by='Abs_NES', ascending=False).head(count)
    print('\nTOP 10 MOST CHANGED PATHWAYS (by Absolute NES):')
    print(df_top[['Term', 'NES', 'FDR q-val']].to_string(index=False))

def run_interactive_explorer():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    h5ad_path = os.path.join(project_dir, 'pik3ca_preprocessed.h5ad')
    try:
        group_a, group_b, suffix = select_comparison()
        pop_key, folder_prefix = select_population()
    except ValueError as error:
        print(f'Error: {error}')
        return
    folder_name = f'{folder_prefix}{suffix}'
    csv_path = os.path.join(script_dir, 'gsea_results', folder_name, 'gseapy.gene_set.prerank.report.csv')
    if not os.path.exists(csv_path):
        print(f'Error: GSEA report not found at {csv_path}')
        return
    df_gsea = pd.read_csv(csv_path)
    show_top_pathways(df_gsea)
    query = input('\nEnter GSEA pathway search term (e.g. Aerobic Respiration): ').strip()
    matched_row = find_matching_pathway(df_gsea, query)
    if matched_row is None:
        print(f"No pathway matched the query: '{query}'")
        return
    term_name = matched_row['Term']
    nes = matched_row['NES']
    lead_genes_str = matched_row['Lead_genes']
    print(f'\nMatched Pathway: {term_name}')
    print(f'Normalized Enrichment Score (NES): {nes:.4f}')
    if pd.isna(lead_genes_str):
        print('No leading edge genes listed for this pathway.')
        return
    genes = str(lead_genes_str).split(';')
    print(f'Found {len(genes)} leading edge genes.')
    print('Loading transcriptomic dataset (this may take a few seconds)...')
    adata = load_dataset(h5ad_path)
    tumor_mask, neuron_mask = compute_population_masks(adata)
    mask = tumor_mask if pop_key == 'tumor' else neuron_mask
    print('Computing differential metrics for driver genes...')
    df_metrics = compute_differential_metrics(adata, mask, group_a, group_b, genes)
    if df_metrics.empty:
        print('No metrics computed (genes might not be present in dataset).')
        return
    df_metrics['Abs_LFC'] = df_metrics['Log2FoldChange'].abs()
    df_sorted = df_metrics.sort_values(by='Abs_LFC', ascending=False).head(20)
    print(f"\nTop 20 Drivers for '{term_name}':")
    print(df_sorted.drop(columns=['Abs_LFC']).to_string(index=False))
if __name__ == '__main__':
    run_interactive_explorer()
