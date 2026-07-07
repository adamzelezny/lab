import pandas as pd
import gseapy as gp

def load_de_results(results_path):
    results_df = pd.read_csv(results_path, index_col=0)
    return results_df


def create_custom_gene_sets():
    gene_sets = {
        "neural_activation": [
            "Fos", "Jun", "Egr1", "Arc", "Bdnf", "Vegfa", "Hif1a", "Cxcl10", "Il6"
        ],
        "cell_cycle_proliferation": [
            "Mki67", "Top2a", "Pcna"
        ],
        "glial_myelin_lineage": [
            "Mog", "Mbp", "Plp1"
        ],
        "astrocyte_markers": [
            "Gfap", "Aqp4"
        ],
        "immune_microglia": [
            "Cd68", "Aif1", "Ccl2", "Cxcl10", "Tnf"
        ],
        "endothelial_vascular": [
            "Pecam1", "Kdr", "Vegfa"
        ]
    }
    return gene_sets


def prepare_ranked_list(results_df):
    ranked_series = results_df["stat"].sort_values(ascending=False)
    return ranked_series


def run_gsea_prerank(ranked_genes, gene_sets):
    pre_res = gp.prerank(
        rnk=ranked_genes,
        gene_sets=gene_sets,
        threads=1,
        min_size=1,
        max_size=100,
        permutation_num=1000,
        outdir=None,
        seed=42
    )
    return pre_res.res2d


def print_gsea_summary(gsea_df):
    summary_cols = ["Term", "ES", "NES", "NOM p-val", "FDR q-val", "Lead_genes"]
    sorted_results = gsea_df[summary_cols].sort_values(by="NES", ascending=False)
    
    print(sorted_results.to_string(index=False))
    print("=========================================================")


def main():
    results_df = load_de_results("differential_expression_results.csv")
    gene_sets = create_custom_gene_sets()
    ranked_genes = prepare_ranked_list(results_df)
    gsea_results_df = run_gsea_prerank(ranked_genes, gene_sets)
    output_filename = "gsea_enrichment_results.csv"
    gsea_results_df.to_csv(output_filename, index=False)
    print(f"\nSaved full GSEA results to '{output_filename}'")
    print_gsea_summary(gsea_results_df)


if __name__ == "__main__":
    main()
