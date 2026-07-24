#Author: Adam Zelezny
#Date: 2026-07-20

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import ranksums

def load_matrix_data(directory_path):
    return sc.read_10x_mtx(directory_path, var_names='gene_symbols')

def process_expression_data(adata):
    sc.pp.normalize_total(adata, target_sum=10000.0)
    sc.pp.log1p(adata)

def compute_gene_metrics(adata, target_genes):
    metrics = {}
    for gene in target_genes:
        expression_values = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()
        expression_values = np.asarray(expression_values)
        non_zero = expression_values > 0
        percentage = float(non_zero.sum() / len(expression_values)) * 100
        mean_val = float(expression_values.mean())
        mean_expressing_val = float(expression_values[non_zero].mean()) if non_zero.sum() > 0 else 0.0
        metrics[gene] = (percentage, mean_val, mean_expressing_val)
    return metrics

def run_significance_tests(ctrl_adata, act_adata, target_genes):
    p_values = {}
    for gene in target_genes:
        ctrl_expr = ctrl_adata[:, gene].X.toarray().flatten() if hasattr(ctrl_adata[:, gene].X, 'toarray') else ctrl_adata[:, gene].X.flatten()
        act_expr = act_adata[:, gene].X.toarray().flatten() if hasattr(act_adata[:, gene].X, 'toarray') else act_adata[:, gene].X.flatten()
        ctrl_expr = np.asarray(ctrl_expr)
        act_expr = np.asarray(act_expr)
        stat, pval = ranksums(ctrl_expr, act_expr)
        p_values[gene] = pval
    return p_values

def save_analysis_results(summary_data, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    summary_data.to_csv(output_path)

def plot_expression_comparison(summary_data, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    summary_data[['Ctrl_Pct%', 'Act_Pct%']].plot(kind='bar', ax=axes[0], color=['#1f77b4', '#ff7f0e'])
    axes[0].set_title('Percentage of Expressing Cells')
    axes[0].set_ylabel('Percentage (%)')
    axes[0].set_xlabel('Genes')
    axes[0].set_xticklabels(summary_data.index, rotation=45)
    summary_data[['Ctrl_Mean', 'Act_Mean']].plot(kind='bar', ax=axes[1], color=['#1f77b4', '#ff7f0e'])
    axes[1].set_title('Mean Expression (All Cells)')
    axes[1].set_ylabel('Mean Expression')
    axes[1].set_xlabel('Genes')
    axes[1].set_xticklabels(summary_data.index, rotation=45)
    summary_data[['Ctrl_Expr_Mean', 'Act_Expr_Mean']].plot(kind='bar', ax=axes[2], color=['#1f77b4', '#ff7f0e'])
    axes[2].set_title('Mean Expression (Expressing Cells)')
    axes[2].set_ylabel('Mean Expression')
    axes[2].set_xlabel('Genes')
    axes[2].set_xticklabels(summary_data.index, rotation=45)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_dot_plot(summary_data, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    genes = summary_data.index.tolist()
    conditions = ['Control', 'Activity']
    x = []
    y = []
    sizes = []
    colors = []
    for idx, gene in enumerate(genes):
        x.append(idx)
        y.append(0)
        sizes.append(summary_data.loc[gene, 'Ctrl_Pct%'])
        colors.append(summary_data.loc[gene, 'Ctrl_Expr_Mean'])
        x.append(idx)
        y.append(1)
        sizes.append(summary_data.loc[gene, 'Act_Pct%'])
        colors.append(summary_data.loc[gene, 'Act_Expr_Mean'])
    fig, ax = plt.subplots(figsize=(14, 5))
    plot_sizes = [s * 10 for s in sizes]
    sc_plot = ax.scatter(x, y, s=plot_sizes, c=colors, cmap='plasma', edgecolors='black', alpha=0.9)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha='right')
    ax.set_yticks([0, 1])
    ax.set_yticklabels(conditions)
    ax.set_ylim(-0.5, 1.5)
    ax.set_xlim(-0.5, len(genes) - 0.5)
    cbar = plt.colorbar(sc_plot, ax=ax, orientation='vertical', shrink=0.75)
    cbar.set_label('Mean Expression (Expressing Cells)')
    handles = []
    labels = []
    for pct in [10, 30, 50, 70]:
        handles.append(ax.scatter([], [], s=pct * 10, c='grey', edgecolors='black', alpha=0.6))
        labels.append(f'{pct}%')
    ax.legend(handles, labels, title='Percent Expressing', loc='upper left', bbox_to_anchor=(1.02, 1.0))
    ax.set_title('TCA Cycle Genes Profile (Dot Plot)', fontsize=14, pad=15)
    ax.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def run_tca_analysis():
    ctrl_dir = 'c:/Users/adamz/Documents/lab/projects/anls/data/control'
    act_dir = 'c:/Users/adamz/Documents/lab/projects/anls/data/activity'
    summary_csv_path = 'c:/Users/adamz/Documents/lab/projects/anls/outputs/tca_summary.csv'
    plot_png_path = 'c:/Users/adamz/Documents/lab/projects/anls/outputs/tca_expression_comparison.png'
    dot_png_path = 'c:/Users/adamz/Documents/lab/projects/anls/outputs/tca_dot_plot.png'
    target_genes = ['Cs', 'Aco2', 'Idh2', 'Idh3a', 'Idh3b', 'Idh3g', 'Ogdh', 'Ogdhl', 'Dlst', 'Dld', 'Suclg1', 'Sucla2', 'Suclg2', 'Sdha', 'Sdhb', 'Sdhc', 'Sdhd', 'Fh1', 'Mdh2']
    ctrl_adata = load_matrix_data(ctrl_dir)
    act_adata = load_matrix_data(act_dir)
    process_expression_data(ctrl_adata)
    process_expression_data(act_adata)
    ctrl_metrics = compute_gene_metrics(ctrl_adata, target_genes)
    act_metrics = compute_gene_metrics(act_adata, target_genes)
    p_values = run_significance_tests(ctrl_adata, act_adata, target_genes)
    rows = []
    for gene in target_genes:
        ctrl_pct, ctrl_mean, ctrl_expr_mean = ctrl_metrics[gene]
        act_pct, act_mean, act_expr_mean = act_metrics[gene]
        pval = p_values[gene]
        rows.append({'Gene': gene, 'Ctrl_Pct%': ctrl_pct, 'Act_Pct%': act_pct, 'Ctrl_Mean': ctrl_mean, 'Act_Mean': act_mean, 'Ctrl_Expr_Mean': ctrl_expr_mean, 'Act_Expr_Mean': act_expr_mean, 'P_Value': pval})
    df = pd.DataFrame(rows).set_index('Gene')
    save_analysis_results(df, summary_csv_path)
    plot_expression_comparison(df, plot_png_path)
    plot_dot_plot(df, dot_png_path)
if __name__ == '__main__':
    run_tca_analysis()
