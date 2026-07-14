import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import MiniBatchKMeans

def load_data(h5ad_path: str) -> sc.AnnData:
    return sc.read_h5ad(h5ad_path)

def preprocess_and_pca(adata: sc.AnnData, n_top_genes: int = 2000, n_comps: int = 50) -> sc.AnnData:
    adata.var['mt'] = adata.var_names.str.lower().str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    n_before = adata.n_obs
    adata = adata[adata.obs['pct_counts_mt'] < 10.0].copy()
    n_after = adata.n_obs
    print(f"Mitochondrial filter (<10%): kept {n_after} of {n_before} cells.")
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=False)
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=n_comps, svd_solver='arpack')
    
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    return adata

def run_kmeans(adata: sc.AnnData, n_clusters: int = 6) -> None:
    X_pca = adata.obsm['X_pca']
    kmeans = MiniBatchKMeans(n_clusters=n_clusters, random_state=42, batch_size=2048)
    clusters = kmeans.fit_predict(X_pca)
    adata.obs['kmeans'] = pd.Series(clusters, index=adata.obs.index).astype(str)
    
    counts = adata.obs['kmeans'].value_counts()
    print("Cluster sizes:")
    for cluster, count in counts.items():
        print(f"  Cluster {cluster}: {count} cells")

def profile_clusters(adata: sc.AnnData, markers: list[str]) -> pd.DataFrame:
    raw_adata = adata.raw.to_adata()
    profile_data = []
    
    for cluster in sorted(adata.obs['kmeans'].unique()):
        cell_mask = adata.obs['kmeans'] == cluster
        cluster_raw = raw_adata[cell_mask]
        
        cluster_stats = {"Cluster": cluster}
        for gene in markers:
            if gene in cluster_raw.var_names:
                gene_idx = cluster_raw.var_names.get_loc(gene)
                if hasattr(cluster_raw.X, "toarray"):
                    expr = cluster_raw.X[:, gene_idx].toarray().flatten()
                else:
                    expr = np.array(cluster_raw.X[:, gene_idx]).flatten()
                
                mean_expr = expr.mean()
                pct_expr = (expr > 0).sum() / len(expr) * 100
                cluster_stats[f"{gene}_mean"] = round(float(mean_expr), 3)
                cluster_stats[f"{gene}_pct"] = round(float(pct_expr), 1)
            else:
                cluster_stats[f"{gene}_mean"] = 0.0
                cluster_stats[f"{gene}_pct"] = 0.0
                
        profile_data.append(cluster_stats)
        
    df_profile = pd.DataFrame(profile_data)
    return df_profile

def plot_pca_results(adata: sc.AnnData, output_dir: str) -> None:
    pca = adata.obsm['X_pca']
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    
    df_plot = pd.DataFrame({
        'PC1': pca[:, 0],
        'PC2': pca[:, 1],
        'Cluster': adata.obs['kmeans'],
        'Variant': adata.obs['variant'],
        'Sample': adata.obs['sample']
    })
    
    raw_adata = adata.raw.to_adata()
    markers_to_plot = ['Gfap', 'Map2']
    for gene in markers_to_plot:
        if gene in raw_adata.var_names:
            gene_idx = raw_adata.var_names.get_loc(gene)
            if hasattr(raw_adata.X, "toarray"):
                expr = raw_adata.X[:, gene_idx].toarray().flatten()
            else:
                expr = np.array(raw_adata.X[:, gene_idx]).flatten()
            df_plot[gene] = expr
            
    if len(df_plot) > 20000:
        df_plot_sub = df_plot.sample(n=20000, random_state=42)
    else:
        df_plot_sub = df_plot
        
    sns.scatterplot(ax=axes[0, 0], data=df_plot_sub, x='PC1', y='PC2', hue='Cluster', 
                    palette='Set1', s=3, alpha=0.6, hue_order=sorted(df_plot['Cluster'].unique()))
    axes[0, 0].set_title('PCA - K-Means Clusters')
    
    sns.scatterplot(ax=axes[0, 1], data=df_plot_sub, x='PC1', y='PC2', hue='Variant', 
                    palette='Dark2', s=3, alpha=0.6)
    axes[0, 1].set_title('PCA - Variants')
    
    if 'Gfap' in df_plot_sub.columns:
        sc1 = axes[1, 0].scatter(df_plot_sub['PC1'], df_plot_sub['PC2'], c=df_plot_sub['Gfap'], 
                                 cmap='viridis', s=2, alpha=0.6)
        fig.colorbar(sc1, ax=axes[1, 0], label='Gfap Expression')
        axes[1, 0].set_title('PCA - Gfap Expression (Glial / Tumor cells)')
        axes[1, 0].set_xlabel('PC1')
        axes[1, 0].set_ylabel('PC2')
        
    if 'Map2' in df_plot_sub.columns:
        sc2 = axes[1, 1].scatter(df_plot_sub['PC1'], df_plot_sub['PC2'], c=df_plot_sub['Map2'], 
                                 cmap='plasma', s=2, alpha=0.6)
        fig.colorbar(sc2, ax=axes[1, 1], label='Map2 Expression')
        axes[1, 1].set_title('PCA - Map2 Expression (Neurons)')
        axes[1, 1].set_xlabel('PC1')
        axes[1, 1].set_ylabel('PC2')
        
    plt.tight_layout()
    plot_path = os.path.join(output_dir, "pca_clustering_results.png")
    plt.savefig(plot_path, dpi=150)
    plt.close()

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    merged_path = os.path.join(project_dir, "pik3ca_merged.h5ad")
    preprocessed_path = os.path.join(project_dir, "pik3ca_preprocessed.h5ad")
    
    if not os.path.exists(merged_path):
        raise FileNotFoundError(f"Merged dataset not found at {merged_path}. Run merge_datasets.py first.")
        
    adata = load_data(merged_path)
    adata = preprocess_and_pca(adata, n_top_genes=2000, n_comps=50)
    run_kmeans(adata, n_clusters=6)
    
    markers = ["Gfap", "Gpc3", "Map2", "C1qa", "Mog", "Flt1"]
    df_profile = profile_clusters(adata, markers)
    print("\nCluster Profiling Table (Mean Expression & Expressing %):")
    print(df_profile.to_string(index=False))
    
    profile_csv_path = os.path.join(script_dir, "cluster_marker_profiles.csv")
    df_profile.to_csv(profile_csv_path, index=False)
    
    cross_tab = pd.crosstab(adata.obs['kmeans'], adata.obs['variant'], normalize='columns') * 100
    print("\nCluster Composition by Variant (% of variant cells in each cluster):")
    print(cross_tab.round(2))
    
    crosstab_csv_path = os.path.join(script_dir, "cluster_variant_composition.csv")
    cross_tab.to_csv(crosstab_csv_path)
    
    plot_pca_results(adata, script_dir)
    
    print(f"\nSaving preprocessed dataset to {preprocessed_path}...")
    adata.write_h5ad(preprocessed_path, compression="gzip")

if __name__ == "__main__":
    main()
