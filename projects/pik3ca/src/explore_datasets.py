#Author: Adam Zelezny
#Date: 2026-07-14

import os
import glob
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def find_samples(data_dir: str) -> list[str]:
    pattern = os.path.join(data_dir, "*_matrix.mtx.gz")
    files = glob.glob(pattern)
    samples = []
    for f in files:
        basename = os.path.basename(f)
        sample = basename.replace("_matrix.mtx.gz", "")
        samples.append(sample)
    return sorted(samples)

def load_sample(data_dir: str, sample_name: str, min_counts: int = 500) -> sc.AnnData:
    adata = sc.read_10x_mtx(data_dir, prefix=f"{sample_name}_", cache=False)
    counts = np.array(adata.X.sum(axis=1)).flatten()
    adata = adata[counts >= min_counts].copy()
    print(f"Loaded {sample_name}: {adata.n_obs} cells")
    return adata

def calculate_qc(adata: sc.AnnData) -> None:
    adata.var['mt'] = adata.var_names.str.lower().str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

def extract_sample_stats(adata: sc.AnnData, sample_name: str) -> dict:
    obs = adata.obs
    active_genes = int((adata.X.sum(axis=0) > 0).sum())
    
    stats = {
        "Sample": sample_name,
        "Cells (>=500 counts)": adata.n_obs,
        "Total Vocabulary Genes": adata.n_vars,
        "Expressed Genes (>=1 count)": active_genes,
        "Mean UMI Counts": float(obs["total_counts"].mean()),
        "Median UMI Counts": float(obs["total_counts"].median()),
        "Std UMI Counts": float(obs["total_counts"].std()),
        "Mean Genes Detected": float(obs["n_genes_by_counts"].mean()),
        "Median Genes Detected": float(obs["n_genes_by_counts"].median()),
        "Std Genes Detected": float(obs["n_genes_by_counts"].std()),
        "Mean Pct Mitochondrial": float(obs["pct_counts_mt"].mean()),
        "Median Pct Mitochondrial": float(obs["pct_counts_mt"].median()),
    }
    return stats

def analyze_key_genes(adata: sc.AnnData, sample_name: str, genes: list[str]) -> list[dict]:
    results = []
    for gene in genes:
        gene_stats = {
            "Sample": sample_name,
            "Gene": gene,
            "Present": False,
            "Pct Expressing Cells": 0.0,
            "Mean Expression (All Cells)": 0.0,
            "Mean Expression (Expressing Cells)": 0.0
        }
        
        if gene in adata.var_names:
            gene_idx = adata.var_names.get_loc(gene)
            if hasattr(adata.X, "toarray"):
                expr = adata.X[:, gene_idx].toarray().flatten()
            else:
                expr = np.array(adata.X[:, gene_idx]).flatten()
                
            gene_stats["Present"] = True
            expressing_mask = expr > 0
            num_expressing = expressing_mask.sum()
            
            if adata.n_obs > 0:
                gene_stats["Pct Expressing Cells"] = float(num_expressing / adata.n_obs * 100)
                gene_stats["Mean Expression (All Cells)"] = float(expr.mean())
            if num_expressing > 0:
                gene_stats["Mean Expression (Expressing Cells)"] = float(expr[expressing_mask].mean())
                
        results.append(gene_stats)
    return results

def plot_qc_metrics(sample_adatas: dict[str, sc.AnnData], output_path: str) -> None:
    obs_dfs = []
    for name, adata in sample_adatas.items():
        df = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
        df["Sample"] = name
        parts = name.split("-")
        df["Condition"] = parts[0]
        obs_dfs.append(df)
        
    combined_df = pd.concat(obs_dfs, ignore_index=True)
    
    fig, axes = plt.subplots(3, 1, figsize=(12, 14), sharex=False)
    
    sns.violinplot(ax=axes[0], data=combined_df, x="Sample", y="total_counts", hue="Condition", palette="Set2")
    axes[0].set_title("Total UMI Counts per Cell")
    axes[0].set_ylabel("UMI Counts")
    axes[0].tick_params(axis='x', rotation=45)
    
    sns.violinplot(ax=axes[1], data=combined_df, x="Sample", y="n_genes_by_counts", hue="Condition", palette="Set2")
    axes[1].set_title("Number of Genes Detected per Cell")
    axes[1].set_ylabel("Genes Detected")
    axes[1].tick_params(axis='x', rotation=45)
    
    sns.violinplot(ax=axes[2], data=combined_df, x="Sample", y="pct_counts_mt", hue="Condition", palette="Set2")
    axes[2].set_title("Mitochondrial Reads Percentage per Cell")
    axes[2].set_ylabel("% Mitochondrial Counts")
    axes[2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()

def df_to_markdown(df: pd.DataFrame, include_index: bool = False) -> str:
    df_copy = df.copy()
    if include_index:
        df_copy = df_copy.reset_index()
    
    cols = df_copy.columns.tolist()
    headers = [str(c) for c in cols]
    
    widths = [len(h) for h in headers]
    rows_data = []
    for _, row in df_copy.iterrows():
        row_cells = [str(row[c]) for c in cols]
        rows_data.append(row_cells)
        for i, val in enumerate(row_cells):
            widths[i] = max(widths[i], len(val))
            
    header_line = "| " + " | ".join(h.ljust(widths[i]) for i, h in enumerate(headers)) + " |"
    separator_line = "| " + " | ".join("-" * widths[i] for i in range(len(cols))) + " |"
    
    table_lines = [header_line, separator_line]
    for row_cells in rows_data:
        row_line = "| " + " | ".join(val.ljust(widths[i]) for i, val in enumerate(row_cells)) + " |"
        table_lines.append(row_line)
        
    return "\n".join(table_lines)

def generate_markdown_report(df_stats: pd.DataFrame, df_genes: pd.DataFrame, report_path: str) -> None:
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# PIK3CA Single-Cell RNA-seq Datasets Summary Report\n\n")
        f.write("This report presents the quality control and expression summary of 10 single-cell RNA-seq datasets ")
        f.write("from the mouse brain glioblastoma model described in Yu et al., Nature 2020. ")
        f.write("The datasets represent replicates of three PIK3CA mutant genotypes: C420R, H1047R, and R88Q.\n\n")
        
        f.write("## Quality Control & Sample Statistics (Cutoff: >= 500 UMI Counts)\n\n")
        f.write("The raw 10x datasets contain empty droplets and GEMs (~6.7M barcodes per sample). ")
        f.write("To isolate real cells, we applied a UMI count threshold of 500 UMI counts. ")
        f.write("The summary statistics of the cells in each sample are shown below:\n\n")
        
        f.write(df_to_markdown(df_stats, include_index=False) + "\n\n")
        
        f.write("## Key Genes of Interest Expression Summary\n\n")
        f.write("Below is the percentage of cells expressing each gene and the mean expression values across cells:\n\n")
        
        df_genes_pivot_pct = df_genes.pivot(index="Gene", columns="Sample", values="Pct Expressing Cells")
        df_genes_pivot_mean = df_genes.pivot(index="Gene", columns="Sample", values="Mean Expression (All Cells)")
        
        f.write("### Percentage of Expressing Cells (%)\n\n")
        f.write(df_to_markdown(df_genes_pivot_pct, include_index=True) + "\n\n")
        
        f.write("### Mean Expression (All Cells)\n\n")
        f.write(df_to_markdown(df_genes_pivot_mean, include_index=True) + "\n\n")
        
        f.write("## Key Findings & Biological Context\n\n")
        f.write("1. **Cell Count Distribution**: The number of cells passing QC (>= 500 UMI counts) ranges across samples, showing variation in cell yield.\n")
        f.write("2. **Driver Expression (Pik3ca)**: Check the expression level of `Pik3ca` to verify presence of PIK3CA mutant driver constructs in IUE-electroporated cells.\n")
        f.write("3. **Tumor Suppressor Check (Pten, Nf1, Trp53)**: These genes are targeted for deletion/knockout to induce tumors. Let's inspect their deletion status by checking if expression is low or absent in cells.\n")
        f.write("4. **Gpc3 (Glypican-3)**: Paper 1 identifies GPC3 as a key driver of gliomagenesis and hyperexcitability. We can examine if Gpc3 expression is higher in C420R/H1047R than in R88Q.\n")
        f.write("5. **Neural vs. Glial markers**: Check the expression of neuronal marker `Map2` and astrocytic markers `Gfap`/`Aldh1l1` to see the composition of the tumor mass and microenvironment.\n")

    print(f"Markdown report generated at {report_path}")

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    data_dir = os.path.join(project_dir, "data")
    
    if not os.path.exists(data_dir):
        raise FileNotFoundError(f"Data directory not found at {data_dir}")
        
    samples = find_samples(data_dir)
    
    genes_of_interest = [
        "Pik3ca",
        "Pten",
        "Nf1",
        "Trp53",
        "Gpc3",
        "Gfap",
        "Aldh1l1",
        "Map2",
        "Dlg4",
        "Syp"
    ]
    
    all_adatas = {}
    stats_list = []
    genes_list = []
    
    for sample in samples:
        adata = load_sample(data_dir, sample, min_counts=500)
        calculate_qc(adata)
        all_adatas[sample] = adata
        
        stats = extract_sample_stats(adata, sample)
        stats_list.append(stats)
        
        gene_stats = analyze_key_genes(adata, sample, genes_of_interest)
        genes_list.extend(gene_stats)
        
    df_stats = pd.DataFrame(stats_list)
    df_genes = pd.DataFrame(genes_list)
    
    stats_csv_path = os.path.join(script_dir, "data_summary.csv")
    genes_csv_path = os.path.join(script_dir, "gene_expression_summary.csv")
    df_stats.to_csv(stats_csv_path, index=False)
    df_genes.to_csv(genes_csv_path, index=False)
    
    plots_path = os.path.join(script_dir, "qc_plots.png")
    plot_qc_metrics(all_adatas, plots_path)
    
    report_path = os.path.join(project_dir, "dataset_summary.md")
    generate_markdown_report(df_stats, df_genes, report_path)

if __name__ == "__main__":
    main()
