#Author: Adam Zelezny
#Date: 2026-07-14

import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, ConfusionMatrixDisplay
import shap

def get_population_masks(adata: sc.AnnData) -> tuple[np.ndarray, np.ndarray]:
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
    print(f'Cells: Tumor={tumor_mask.sum()}, Neurons={neuron_mask.sum()}')
    return (tumor_mask, neuron_mask)

def train_and_explain(adata: sc.AnnData, mask: np.ndarray, pop_name: str, out_dir: str) -> dict[str, pd.DataFrame]:
    print(f'\n--- Processing {pop_name} Classification ---')
    adata_sub = adata[mask].copy()
    raw_adata = adata.raw.to_adata()[mask].copy()
    hvg_mask = adata_sub.var['highly_variable']
    hvg_names = adata_sub.var_names[hvg_mask].tolist()
    X = raw_adata[:, hvg_names].X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    else:
        X = np.array(X)
    X_df = pd.DataFrame(X, columns=hvg_names, index=adata_sub.obs_names)
    y = adata_sub.obs['variant'].values
    print(f'Class distribution: {pd.Series(y).value_counts().to_dict()}')
    X_train, X_test, y_train, y_test = train_test_split(X_df, y, test_size=0.2, stratify=y, random_state=42)
    rf = RandomForestClassifier(n_estimators=100, max_depth=10, class_weight='balanced', random_state=42, n_jobs=-1)
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)
    print(f'\nClassification Report for {pop_name}:')
    print(classification_report(y_test, y_pred))
    fig, ax = plt.subplots(figsize=(6, 5))
    cm = confusion_matrix(y_test, y_pred, labels=rf.classes_)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=rf.classes_)
    disp.plot(cmap='Blues', ax=ax)
    plt.title(f'{pop_name} Classification CM')
    plt.tight_layout()
    cm_path = os.path.join(out_dir, f'{pop_name.lower()}_confusion_matrix.png')
    plt.savefig(cm_path, dpi=150)
    plt.close()
    print('Calculating SHAP values...')
    max_shap_samples = 200
    if len(X_test) > max_shap_samples:
        X_test_sub = X_test.sample(n=max_shap_samples, random_state=42)
    else:
        X_test_sub = X_test
    explainer = shap.TreeExplainer(rf)
    shap_values = explainer(X_test_sub)
    top_drivers = {}
    for k, class_name in enumerate(rf.classes_):
        mean_abs_shap = np.abs(shap_values.values[:, :, k]).mean(axis=0)
        df_shap = pd.DataFrame({'Gene': X_df.columns, 'Mean_Abs_SHAP': mean_abs_shap})
        df_shap = df_shap.sort_values(by='Mean_Abs_SHAP', ascending=False).head(15)
        top_drivers[class_name] = df_shap
        plt.figure(figsize=(8, 6))
        shap.plots.bar(shap_values[:, :, k], max_display=15, show=False)
        plt.title(f'{pop_name} - Top Drivers for {class_name}')
        plt.tight_layout()
        bar_path = os.path.join(out_dir, f'{pop_name.lower()}_shap_bar_{class_name.lower()}.png')
        plt.savefig(bar_path, dpi=150)
        plt.close()
        plt.figure(figsize=(8, 6))
        shap.plots.beeswarm(shap_values[:, :, k], max_display=15, show=False)
        plt.title(f'{pop_name} - Beeswarm for {class_name}')
        plt.tight_layout()
        beeswarm_path = os.path.join(out_dir, f'{pop_name.lower()}_shap_beeswarm_{class_name.lower()}.png')
        plt.savefig(beeswarm_path, dpi=150)
        plt.close()
    return top_drivers

def df_to_markdown(df: pd.DataFrame) -> str:
    cols = df.columns.tolist()
    headers = [str(c) for c in cols]
    widths = [len(h) for h in headers]
    rows_data = []
    for _, row in df.iterrows():
        row_cells = [str(row[c]) for c in cols]
        rows_data.append(row_cells)
        for i, val in enumerate(row_cells):
            widths[i] = max(widths[i], len(val))
    header_line = '| ' + ' | '.join((h.ljust(widths[i]) for i, h in enumerate(headers))) + ' |'
    separator_line = '| ' + ' | '.join(('-' * widths[i] for i in range(len(cols)))) + ' |'
    table_lines = [header_line, separator_line]
    for row_cells in rows_data:
        row_line = '| ' + ' | '.join((val.ljust(widths[i]) for i, val in enumerate(row_cells))) + ' |'
        table_lines.append(row_line)
    return '\n'.join(table_lines)

def write_summary_report(tumor_drivers: dict[str, pd.DataFrame], neuron_drivers: dict[str, pd.DataFrame], report_path: str) -> None:
    print(f'Writing SHAP report to {report_path}...')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('# Single-Cell Genotype Classification & SHAP Analysis Report\n\n')
        f.write('This report presents the machine learning classification and feature importance (SHAP) results ')
        f.write('for individual cells. We trained **Random Forest Classifiers** to distinguish between `R88Q`, `C420R`, ')
        f.write('and `H1047R` genotypes, capturing non-linear interactions at the single-cell level.\n\n')
        f.write('## 1. Malignant Tumor Cells SHAP Analysis\n\n')
        for genotype, df in tumor_drivers.items():
            f.write(f'### Top 15 Driver Genes for Genotype: **{genotype}**\n\n')
            f.write(df_to_markdown(df) + '\n\n')
        f.write('## 2. Microenvironmental Neurons SHAP Analysis\n\n')
        for genotype, df in neuron_drivers.items():
            f.write(f'### Top 15 Driver Genes for Genotype: **{genotype}**\n\n')
            f.write(df_to_markdown(df) + '\n\n')

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    preprocessed_path = os.path.join(project_dir, 'pik3ca_preprocessed.h5ad')
    if not os.path.exists(preprocessed_path):
        raise FileNotFoundError(f'Preprocessed file not found at {preprocessed_path}.')
    adata = sc.read_h5ad(preprocessed_path)
    tumor_mask, neuron_mask = get_population_masks(adata)
    tumor_drivers = train_and_explain(adata, tumor_mask, 'Tumor', script_dir)
    neuron_drivers = train_and_explain(adata, neuron_mask, 'Neuron', script_dir)
    report_path = os.path.join(project_dir, 'shap_summary.md')
    write_summary_report(tumor_drivers, neuron_drivers, report_path)
    print('SHAP pipeline execution complete!')
if __name__ == '__main__':
    main()
