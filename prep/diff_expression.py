import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def load_and_align_data(counts_path, metadata_path):
    counts_df = pd.read_csv(counts_path, index_col=0)
    metadata_df = pd.read_csv(metadata_path)    
    metadata_df.set_index("cell_id", inplace=True)    
    metadata_df = metadata_df.loc[counts_df.index]
    initial_genes = counts_df.shape[1]
    real_genes = [col for col in counts_df.columns if not col.startswith("Gene_")]
    counts_df = counts_df[real_genes]    
    return counts_df, metadata_df


def run_deseq2(counts_df, metadata_df, design_factor="condition"):
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        design=f"~{design_factor}",
        size_factors_fit_type="poscounts",
        refit_cooks=True,
        n_cpus=1
    )
    dds.deseq2()
    return dds


def get_deseq_results(dds, contrast):
    stat_res = DeseqStats(
        dds,
        contrast=contrast,
        n_cpus=1
    )
    stat_res.summary()
    return stat_res.results_df

def print_results_summary(results_df, padj_threshold=0.05, lfc_threshold=1.0):
    #filter for statistically significant genes
    significant_genes = results_df[results_df["padj"] < padj_threshold]
    
    print(f"\n=== Top Upregulated Genes (padj < {padj_threshold}, log2FC > {lfc_threshold}) ===")
    upregulated = significant_genes[significant_genes["log2FoldChange"] > lfc_threshold].sort_values(
        by="padj"
    )
    print(upregulated[["log2FoldChange", "pvalue", "padj"]].head(15))
    
    print(f"\n=== Top Downregulated Genes (padj < {padj_threshold}, log2FC < -{lfc_threshold}) ===")
    downregulated = significant_genes[significant_genes["log2FoldChange"] < -lfc_threshold].sort_values(
        by="padj"
    )
    print(downregulated[["log2FoldChange", "pvalue", "padj"]].head(15))


def main():
    # 1.load and align data
    counts_df, metadata_df = load_and_align_data(
        "synthetic_glioblastoma_counts_500x500.csv",
        "synthetic_glioblastoma_metadata.csv"
    )
    
    # 2.fit the DESeq2 model
    dds = run_deseq2(counts_df, metadata_df, design_factor="condition")
    
    # 3.perform statistical tests (activation vs pre_activation)
    results_df = get_deseq_results(
        dds, 
        contrast=["condition", "activation", "pre_activation"]
    )
    
    # 4.save results to disk
    output_filename = "differential_expression_results.csv"
    results_df.to_csv(output_filename)
    print(f"\nSaved full results to '{output_filename}'")
    
    # 5.summarize the differential expression
    print_results_summary(results_df, padj_threshold=0.05, lfc_threshold=1.0)


if __name__ == "__main__":
    main()
