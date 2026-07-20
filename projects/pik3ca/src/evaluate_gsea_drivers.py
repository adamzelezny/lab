#Author: Adam Zelezny
#Date: 2026-07-17

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.stats


def load_data_matrix(h5ad_path):
    return sc.read_h5ad(h5ad_path)


def extract_population_masks(adata):
    raw_adata = adata.raw.to_adata()
    gfap_idx = raw_adata.var_names.get_loc("Gfap")
    map2_idx = raw_adata.var_names.get_loc("Map2")
    gpc3_idx = raw_adata.var_names.get_loc("Gpc3")
    
    if hasattr(raw_adata.X, "toarray"):
        gfap_expr = raw_adata.X[:, gfap_idx].toarray().flatten()
        map2_expr = raw_adata.X[:, map2_idx].toarray().flatten()
        gpc3_expr = raw_adata.X[:, gpc3_idx].toarray().flatten()
    else:
        gfap_expr = np.array(raw_adata.X[:, gfap_idx]).flatten()
        map2_expr = np.array(raw_adata.X[:, map2_idx]).flatten()
        gpc3_expr = np.array(raw_adata.X[:, gpc3_idx]).flatten()
        
    tumor_mask = (adata.obs["kmeans"] == "2") & ((gfap_expr > 0) | (gpc3_expr > 0)) & (map2_expr == 0)
    neuron_mask = (map2_expr > 0) & (gfap_expr == 0)
    return tumor_mask, neuron_mask


def parse_gsea_leading_edge(csv_path, top_n_pathways=3, top_n_genes=10):
    if not os.path.exists(csv_path):
        return []
        
    df = pd.read_csv(csv_path)
    
    upregulated_df = df[df["NES"] > 0].sort_values(by="NES", ascending=False)
    downregulated_df = df[df["NES"] < 0].sort_values(by="NES", ascending=True)
    
    pathways = []
    
    for _, row in upregulated_df.head(top_n_pathways).iterrows():
        term = row["Term"]
        nes = row["NES"]
        if pd.isna(row["Lead_genes"]):
            continue
        genes = str(row["Lead_genes"]).split(";")[:top_n_genes]
        pathways.append((term, nes, "Upregulated", genes))
        
    for _, row in downregulated_df.head(top_n_pathways).iterrows():
        term = row["Term"]
        nes = row["NES"]
        if pd.isna(row["Lead_genes"]):
            continue
        genes = str(row["Lead_genes"]).split(";")[:top_n_genes]
        pathways.append((term, nes, "Downregulated", genes))
        
    return pathways


def compute_differential_metrics(adata, mask, group_a, group_b, genes):
    raw_adata = adata.raw.to_adata()[mask]
    
    results = []
    for gene in genes:
        if gene not in raw_adata.var_names:
            continue
            
        gene_idx = raw_adata.var_names.get_loc(gene)
        
        cells_a = (raw_adata.obs["variant"] == group_a).values
        cells_b = (raw_adata.obs["variant"] == group_b).values
        
        if hasattr(raw_adata.X, "toarray"):
            expr_a = raw_adata.X[cells_a, gene_idx].toarray().flatten()
            expr_b = raw_adata.X[cells_b, gene_idx].toarray().flatten()
        else:
            expr_a = np.array(raw_adata.X[cells_a, gene_idx]).flatten()
            expr_b = np.array(raw_adata.X[cells_b, gene_idx]).flatten()
            
        mean_a = float(np.mean(expr_a))
        mean_b = float(np.mean(expr_b))
        
        pct_a = float(np.mean(expr_a > 0) * 100)
        pct_b = float(np.mean(expr_b > 0) * 100)
        
        log2_fc = float(np.log2(mean_a + 1e-4) - np.log2(mean_b + 1e-4))
        
        if len(expr_a) > 0 and len(expr_b) > 0:
            stat, p_val = scipy.stats.mannwhitneyu(expr_a, expr_b, alternative="two-sided")
        else:
            p_val = 1.0
            
        results.append({
            "Gene": gene,
            f"Mean_{group_a}": mean_a,
            f"Mean_{group_b}": mean_b,
            f"Pct_{group_a}": pct_a,
            f"Pct_{group_b}": pct_b,
            "Log2FoldChange": log2_fc,
            "PValue": p_val
        })
        
    return pd.DataFrame(results)


def run_driver_evaluation():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    h5ad_path = os.path.join(project_dir, "pik3ca_preprocessed.h5ad")
    
    adata = load_data_matrix(h5ad_path)
    tumor_mask, neuron_mask = extract_population_masks(adata)
    
    comparisons = [
        ("tumor", "H1047R", "R88Q", "tumor_H1047R_vs_R88Q_control", tumor_mask),
        ("tumor", "C420R", "R88Q", "tumor_C420R_vs_R88Q_control", tumor_mask),
        ("tumor", "H1047R", "C420R", "tumor_H1047R_vs_C420R_direct", tumor_mask),
        ("neuron", "H1047R", "R88Q", "neuron_H1047R_vs_R88Q_control", neuron_mask),
        ("neuron", "C420R", "R88Q", "neuron_C420R_vs_R88Q_control", neuron_mask),
        ("neuron", "H1047R", "C420R", "neuron_H1047R_vs_C420R_direct", neuron_mask)
    ]
    
    for pop_name, group_a, group_b, folder, mask in comparisons:
        csv_path = os.path.join(script_dir, "gsea_results", folder, "gseapy.gene_set.prerank.report.csv")
        if not os.path.exists(csv_path):
            print(f"Skipping {folder}: report not found.")
            continue
            
        print(f"\nEvaluating drivers for: {folder}")
        pathway_genes = parse_gsea_leading_edge(csv_path, top_n_pathways=3, top_n_genes=10)
        
        all_dfs = []
        for term, nes, direction, genes in pathway_genes:
            df = compute_differential_metrics(adata, mask, group_a, group_b, genes)
            if not df.empty:
                df["Pathway"] = term
                df["NES"] = nes
                df["Direction"] = direction
                all_dfs.append(df)
                
        if all_dfs:
            final_df = pd.concat(all_dfs, ignore_index=True)
            output_csv_path = os.path.join(script_dir, "gsea_results", folder, "leading_edge_drivers_raw.csv")
            final_df.to_csv(output_csv_path, index=False)
            print(f"Saved results to {output_csv_path}")
            
            final_df["Abs_LFC"] = final_df["Log2FoldChange"].abs()
            top_drivers = final_df.sort_values(by="Abs_LFC", ascending=False).head(5)
            print("Top 5 Driver Genes by Log2FoldChange:")
            for _, row in top_drivers.iterrows():
                print(f"  {row['Gene']} ({row['Pathway']}): LFC={row['Log2FoldChange']:.4f}, p={row['PValue']:.2e}")


if __name__ == "__main__":
    run_driver_evaluation()
