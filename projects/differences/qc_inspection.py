#Author: Adam Zelezny
#Date: 2026-07-08

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
DATA_FILE = Path('projects/differences/data/combined_labeled_data.h5ad')
OUT_DIR = Path('projects/differences/outputs/qc')
OUT_DIR.mkdir(parents=True, exist_ok=True)

def print_section(title: str):
    print('\n' + '=' * 70)
    print(title)
    print('=' * 70)

def summarize_numeric(series: pd.Series, label: str):
    print(f'\n--- {label} ---')
    print(f'Mean   : {series.mean():.2f}')
    print(f'Median : {series.median():.2f}')
    print(f'Std    : {series.std():.2f}')
    print(f'Min    : {series.min():.2f}')
    print(f'1%     : {series.quantile(0.01):.2f}')
    print(f'5%     : {series.quantile(0.05):.2f}')
    print(f'95%    : {series.quantile(0.95):.2f}')
    print(f'99%    : {series.quantile(0.99):.2f}')
    print(f'Max    : {series.max():.2f}')

def save_histogram(adata, column, filename, title, xlabel, bins=80):
    plt.figure(figsize=(8, 5))
    for condition in adata.obs['condition'].unique():
        values = adata.obs.loc[adata.obs['condition'] == condition, column]
        plt.hist(values, bins=bins, alpha=0.5, label=condition)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('Number of cells')
    plt.legend()
    plt.tight_layout()
    out_path = OUT_DIR / filename
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f'Saved: {out_path}')

def save_scatter(adata, x_col, y_col, filename, title, xlabel, ylabel, max_points=20000):
    obs = adata.obs.copy()
    if len(obs) > max_points:
        obs = obs.sample(max_points, random_state=42)
    plt.figure(figsize=(7, 5))
    for condition in obs['condition'].unique():
        subset = obs[obs['condition'] == condition]
        plt.scatter(subset[x_col], subset[y_col], s=4, alpha=0.35, label=condition)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(markerscale=3)
    plt.tight_layout()
    out_path = OUT_DIR / filename
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f'Saved: {out_path}')

def save_violin(adata, columns, filename):
    sc.pl.violin(adata, keys=columns, groupby='condition', jitter=0.4, multi_panel=True, show=False)
    out_path = OUT_DIR / filename
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')
print_section('LOAD COMBINED DATASET')
print(f'Loading: {DATA_FILE}')
adata = sc.read_h5ad(DATA_FILE)
print(adata)
print('\nCell metadata columns:')
print(list(adata.obs.columns))
print('\nGene metadata columns:')
print(list(adata.var.columns))
print_section('IDENTIFY QC GENE GROUPS')
adata.var['mt'] = adata.var_names.str.startswith('mt-')
adata.var['ribo'] = adata.var_names.str.startswith(('Rps', 'Rpl'))
print(f"Mitochondrial genes detected: {int(adata.var['mt'].sum()):,}")
print(f"Ribosomal genes detected    : {int(adata.var['ribo'].sum()):,}")
print('\nFirst mitochondrial genes:')
print(adata.var_names[adata.var['mt']][:10].tolist())
print('\nFirst ribosomal genes:')
print(adata.var_names[adata.var['ribo']][:10].tolist())
print_section('COMPUTE QC METRICS')
print('Computing total counts, detected genes, mitochondrial %, and ribosomal %...')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
print('QC metrics added to adata.obs.')
print('\nNew QC columns:')
qc_columns = ['total_counts', 'n_genes_by_counts', 'total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'pct_counts_ribo']
print(qc_columns)
print_section('OVERALL QC SUMMARY')
print(f'Total cells: {adata.n_obs:,}')
print(f'Total genes: {adata.n_vars:,}')
summarize_numeric(adata.obs['total_counts'], 'Total RNA counts per cell')
summarize_numeric(adata.obs['n_genes_by_counts'], 'Genes detected per cell')
summarize_numeric(adata.obs['pct_counts_mt'], 'Mitochondrial percent per cell')
summarize_numeric(adata.obs['pct_counts_ribo'], 'Ribosomal percent per cell')
print_section('QC SUMMARY BY CONDITION')
summary = adata.obs.groupby('condition').agg(cells=('condition', 'count'), mean_counts=('total_counts', 'mean'), median_counts=('total_counts', 'median'), mean_genes=('n_genes_by_counts', 'mean'), median_genes=('n_genes_by_counts', 'median'), mean_pct_mt=('pct_counts_mt', 'mean'), median_pct_mt=('pct_counts_mt', 'median'), mean_pct_ribo=('pct_counts_ribo', 'mean'), median_pct_ribo=('pct_counts_ribo', 'median'))
print(summary.round(2))
summary_file = OUT_DIR / 'qc_summary_by_condition.csv'
summary.to_csv(summary_file)
print(f'\nSaved summary table: {summary_file}')
print_section('POSSIBLE QC THRESHOLD CHECKS')
candidate_thresholds = {'low_gene_cells__n_genes_lt_200': adata.obs['n_genes_by_counts'] < 200, 'low_gene_cells__n_genes_lt_500': adata.obs['n_genes_by_counts'] < 500, 'high_gene_cells__n_genes_gt_6000': adata.obs['n_genes_by_counts'] > 6000, 'high_gene_cells__n_genes_gt_7000': adata.obs['n_genes_by_counts'] > 7000, 'high_mito_cells__pct_mt_gt_10': adata.obs['pct_counts_mt'] > 10, 'high_mito_cells__pct_mt_gt_15': adata.obs['pct_counts_mt'] > 15, 'high_mito_cells__pct_mt_gt_20': adata.obs['pct_counts_mt'] > 20, 'low_count_cells__counts_lt_1000': adata.obs['total_counts'] < 1000, 'high_count_cells__counts_gt_30000': adata.obs['total_counts'] > 30000}
threshold_rows = []
for name, mask in candidate_thresholds.items():
    total_failed = int(mask.sum())
    pct_failed = 100 * total_failed / adata.n_obs
    row = {'criterion': name, 'cells_failed': total_failed, 'percent_failed': pct_failed}
    for condition in adata.obs['condition'].unique():
        condition_mask = adata.obs['condition'] == condition
        failed_in_condition = int((mask & condition_mask).sum())
        condition_total = int(condition_mask.sum())
        row[f'{condition}_failed'] = failed_in_condition
        row[f'{condition}_percent_failed'] = 100 * failed_in_condition / condition_total
    threshold_rows.append(row)
threshold_df = pd.DataFrame(threshold_rows)
print(threshold_df.round(2).to_string(index=False))
threshold_file = OUT_DIR / 'candidate_threshold_failures.csv'
threshold_df.to_csv(threshold_file, index=False)
print(f'\nSaved threshold report: {threshold_file}')
print_section('SAVE PER-CELL QC METRICS')
qc_table = adata.obs[['condition', 'mouse_id', 'timepoint', 'source_folder', 'dataset', 'total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribo']].copy()
qc_table_file = OUT_DIR / 'per_cell_qc_metrics.csv'
qc_table.to_csv(qc_table_file)
print(f'Saved per-cell QC table: {qc_table_file}')
print_section('SAVE QC PLOTS')
save_histogram(adata, column='total_counts', filename='hist_total_counts_by_condition.png', title='Total RNA counts per cell', xlabel='Total counts')
save_histogram(adata, column='n_genes_by_counts', filename='hist_genes_detected_by_condition.png', title='Genes detected per cell', xlabel='Number of genes detected')
save_histogram(adata, column='pct_counts_mt', filename='hist_mito_percent_by_condition.png', title='Mitochondrial percentage per cell', xlabel='Percent mitochondrial counts')
save_histogram(adata, column='pct_counts_ribo', filename='hist_ribo_percent_by_condition.png', title='Ribosomal percentage per cell', xlabel='Percent ribosomal counts')
save_scatter(adata, x_col='total_counts', y_col='n_genes_by_counts', filename='scatter_counts_vs_genes.png', title='Total counts vs genes detected', xlabel='Total counts', ylabel='Genes detected')
save_scatter(adata, x_col='total_counts', y_col='pct_counts_mt', filename='scatter_counts_vs_mito.png', title='Total counts vs mitochondrial percentage', xlabel='Total counts', ylabel='Percent mitochondrial counts')
save_violin(adata, columns=['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribo'], filename='violin_qc_metrics_by_condition.png')
print_section('SAVE QC-ANNOTATED DATASET')
qc_h5ad = OUT_DIR / 'combined_labeled_with_qc_metrics.h5ad'
adata.write_h5ad(qc_h5ad)
print(f'Saved QC-annotated dataset: {qc_h5ad}')
print('\nQC inspection complete.')
print('No cells were filtered in this script.')
print('Review the plots and threshold report before choosing final filtering thresholds.')
