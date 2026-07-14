import os
import scanpy as sc
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp

def get_population_masks(adata: sc.AnnData) -> tuple[np.ndarray, np.ndarray]:
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
        
    tumor_mask = (adata.obs['kmeans'] == '2') & ((gfap_expr > 0) | (gpc3_expr > 0)) & (map2_expr == 0)
    neuron_mask = (map2_expr > 0) & (gfap_expr == 0)
    
    print(f"Purified Tumor Cells: {tumor_mask.sum()}")
    print(f"Purified Neurons: {neuron_mask.sum()}")
    return tumor_mask, neuron_mask

def make_pseudobulk(adata_raw: sc.AnnData, cell_mask: np.ndarray, sample_col: str = 'sample') -> tuple[pd.DataFrame, pd.DataFrame]:
    adata_sub = adata_raw[cell_mask]
    X = adata_sub.X
    samples = adata_sub.obs[sample_col]
    unique_samples = sorted(samples.unique())
    
    pseudobulk_matrix = np.zeros((len(unique_samples), adata_raw.n_vars))
    for i, s in enumerate(unique_samples):
        s_mask = (samples == s).values
        s_counts = X[s_mask].sum(axis=0)
        pseudobulk_matrix[i, :] = np.array(s_counts).flatten()
        
    df_counts = pd.DataFrame(
        pseudobulk_matrix,
        index=unique_samples,
        columns=adata_raw.var_names
    )
    df_counts = df_counts.astype(int)
    
    genes_to_keep = df_counts.sum(axis=0) >= 10
    df_counts = df_counts.loc[:, genes_to_keep]
    
    df_metadata = pd.DataFrame(index=unique_samples)
    df_metadata['variant'] = [s.split('-')[0] for s in unique_samples]
    df_metadata['replicate'] = [s.split('-')[1] for s in unique_samples]
    
    return df_counts, df_metadata

def run_deseq2(counts_df: pd.DataFrame, metadata_df: pd.DataFrame, out_dir: str, prefix: str) -> dict[str, pd.DataFrame]:
    print(f"Running PyDESeq2 for {prefix}...")
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        design_factors="variant",
        refit_cooks=True
    )
    dds.deseq2()
    
    results = {}
    contrasts = [
        ("H1047R", "R88Q", "H1047R_vs_R88Q_control"),
        ("C420R", "R88Q", "C420R_vs_R88Q_control"),
        ("H1047R", "C420R", "H1047R_vs_C420R_direct")
    ]
    
    for num, denom, label in contrasts:
        stat_res = DeseqStats(dds, contrast=["variant", num, denom])
        stat_res.summary()
        res_df = stat_res.results_df
        
        csv_name = f"DE_{prefix}_{label}.csv"
        csv_path = os.path.join(out_dir, csv_name)
        res_df.to_csv(csv_path)
        results[label] = res_df
        
    return results

def run_gsea(de_results: dict[str, pd.DataFrame], output_dir: str, prefix: str) -> dict[str, pd.DataFrame]:
    print(f"Running GSEA for {prefix}...")
    gsea_summaries = {}
    
    for label, res_df in de_results.items():
        res_clean = res_df.dropna(subset=['stat'])
        ranking = res_clean['stat'].sort_values(ascending=False)
        
        contrast_out_dir = os.path.join(output_dir, "gsea_results", f"{prefix}_{label}")
        os.makedirs(contrast_out_dir, exist_ok=True)
        
        try:
            pre_res = gp.prerank(
                rnk=ranking,
                gene_sets='GO_Biological_Process_2023',
                outdir=contrast_out_dir,
                permutation_num=100,
                min_size=15,
                max_size=500,
                no_plot=True
            )
            gsea_df = pre_res.res2d
            gsea_summaries[label] = gsea_df.sort_values(by='NES', key=abs, ascending=False).head(15)
        except Exception as e:
            print(f"  [Warning] GSEA failed for {label}: {e}")
            gsea_summaries[label] = pd.DataFrame()
            
    return gsea_summaries

def df_to_markdown_custom(df: pd.DataFrame, include_index: bool = False) -> str:
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

def generate_summary_report(
    tumor_de: dict[str, pd.DataFrame], 
    neuron_de: dict[str, pd.DataFrame],
    tumor_gsea: dict[str, pd.DataFrame],
    neuron_gsea: dict[str, pd.DataFrame],
    report_path: str
) -> None:
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# PIK3CA Variant Differential Expression & Pathway Analysis Report\n\n")
        f.write("This report presents the differential expression (PyDESeq2) and pathway enrichment (GSEA Pre-ranked) results ")
        f.write("for purified **Tumor Cells** and **Neurons** across the PIK3CA mutant genotypes. ")
        f.write("The variant **R88Q** is treated as the control baseline group (no activity) to highlight the activated properties ")
        f.write("of **H1047R** and **C420R**.\n\n")
        
        f.write("## 1. Malignant Tumor Cells Analysis\n\n")
        f.write("Tumor cells were defined as cells from Cluster 5 expressing glial/tumor markers (`Gfap` or `Gpc3`) but lacking neuronal marker (`Map2`). ")
        f.write("Comparing these cells isolates the signals that each tumor type sends to the brain microenvironment.\n\n")
        
        for label in ["H1047R_vs_R88Q_control", "C420R_vs_R88Q_control", "H1047R_vs_C420R_direct"]:
            f.write(f"### Contrast: {label.replace('_', ' ')}\n\n")
            
            de_df = tumor_de[label]
            up = de_df[de_df['padj'] < 0.05].sort_values(by='log2FoldChange', ascending=False).head(10)[['baseMean', 'log2FoldChange', 'padj']]
            down = de_df[de_df['padj'] < 0.05].sort_values(by='log2FoldChange', ascending=True).head(10)[['baseMean', 'log2FoldChange', 'padj']]
            
            f.write("#### Top 10 Upregulated Genes (Significance padj < 0.05):\n\n")
            if not up.empty:
                f.write(df_to_markdown_custom(up, include_index=True) + "\n\n")
            else:
                f.write("*No significant genes found.*\n\n")
                
            f.write("#### Top 10 Downregulated Genes (Significance padj < 0.05):\n\n")
            if not down.empty:
                f.write(df_to_markdown_custom(down, include_index=True) + "\n\n")
            else:
                f.write("*No significant genes found.*\n\n")
                
            gsea_df = tumor_gsea[label]
            if not gsea_df.empty:
                f.write("#### Top Enriched Biological Process Pathways (NES sorted):\n\n")
                f.write(df_to_markdown_custom(gsea_df[['Term', 'NES', 'NOM p-val', 'FDR q-val']], include_index=False) + "\n\n")
                
        f.write("## 2. Microenvironmental Neurons Analysis\n\n")
        f.write("Neurons were defined as cells expressing neuronal markers (`Map2` or `Dlg4`) but lacking glial marker (`Gfap`). ")
        f.write("Comparing these cells reveals the transcription changes occurring in host neurons in response to the surrounding tumor variants.\n\n")
        
        for label in ["H1047R_vs_R88Q_control", "C420R_vs_R88Q_control", "H1047R_vs_C420R_direct"]:
            f.write(f"### Contrast: {label.replace('_', ' ')}\n\n")
            
            de_df = neuron_de[label]
            up = de_df[de_df['padj'] < 0.05].sort_values(by='log2FoldChange', ascending=False).head(10)[['baseMean', 'log2FoldChange', 'padj']]
            down = de_df[de_df['padj'] < 0.05].sort_values(by='log2FoldChange', ascending=True).head(10)[['baseMean', 'log2FoldChange', 'padj']]
            
            f.write("#### Top 10 Upregulated Genes (Significance padj < 0.05):\n\n")
            if not up.empty:
                f.write(df_to_markdown_custom(up, include_index=True) + "\n\n")
            else:
                f.write("*No significant genes found.*\n\n")
                
            f.write("#### Top 10 Downregulated Genes (Significance padj < 0.05):\n\n")
            if not down.empty:
                f.write(df_to_markdown_custom(down, include_index=True) + "\n\n")
            else:
                f.write("*No significant genes found.*\n\n")
                
            gsea_df = neuron_gsea[label]
            if not gsea_df.empty:
                f.write("#### Top Enriched Biological Process Pathways (NES sorted):\n\n")
                f.write(df_to_markdown_custom(gsea_df[['Term', 'NES', 'NOM p-val', 'FDR q-val']], include_index=False) + "\n\n")

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    
    preprocessed_path = os.path.join(project_dir, "pik3ca_preprocessed.h5ad")
    merged_path = os.path.join(project_dir, "pik3ca_merged.h5ad")
    
    if not os.path.exists(preprocessed_path):
        raise FileNotFoundError(f"Preprocessed file not found at {preprocessed_path}.")
    if not os.path.exists(merged_path):
        raise FileNotFoundError(f"Merged raw file not found at {merged_path}.")
        
    adata_prep = sc.read_h5ad(preprocessed_path)
    adata_raw = sc.read_h5ad(merged_path)
    
    adata_raw = adata_raw[adata_prep.obs_names].copy()
    adata_raw.obs['kmeans'] = adata_prep.obs['kmeans']
    
    tumor_mask, neuron_mask = get_population_masks(adata_prep)
    
    tumor_counts, tumor_metadata = make_pseudobulk(adata_raw, tumor_mask)
    neuron_counts, neuron_metadata = make_pseudobulk(adata_raw, neuron_mask)
    
    tumor_de_results = run_deseq2(tumor_counts, tumor_metadata, script_dir, "tumor")
    neuron_de_results = run_deseq2(neuron_counts, neuron_metadata, script_dir, "neuron")
    
    tumor_gsea = run_gsea(tumor_de_results, script_dir, "tumor")
    neuron_gsea = run_gsea(neuron_de_results, script_dir, "neuron")
    
    report_path = os.path.join(project_dir, "de_gsea_summary.md")
    generate_summary_report(tumor_de_results, neuron_de_results, tumor_gsea, neuron_gsea, report_path)

if __name__ == "__main__":
    main()
