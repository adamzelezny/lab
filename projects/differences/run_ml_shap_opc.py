#Author: Adam Zelezny
#Date: 2026-07-09

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
import shap
import matplotlib.pyplot as plt
INPUT_FILE = Path('projects/differences/outputs/annotated/combined_annotated.h5ad')
OUT_DIR = Path('projects/differences/outputs/ml_shap')
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_CSV = OUT_DIR / 'opc_activity_shap_candidates.csv'
PLOT_FILE = OUT_DIR / 'opc_activity_shap_summary.png'
print('Loading annotated dataset...')
adata = sc.read_h5ad(INPUT_FILE)
print(adata)
print(f'\nFiltering dataset to OPC-like cells (Cluster 12)...')
adata_opc = adata[adata.obs['cell_type'] == 'OPC-like cells'].copy()
print(adata_opc)
print('\nCell counts by condition in OPC-like subset:')
print(adata_opc.obs['condition'].value_counts())
highly_variable_genes_mask = adata_opc.var['highly_variable']
feature_gene_names = adata_opc.var_names[highly_variable_genes_mask].tolist()
print(f'\nExtracting expression values for {len(feature_gene_names)} highly variable genes...')
X = adata_opc[:, highly_variable_genes_mask].X.toarray()
y = (adata_opc.obs['condition'] == 'activity').astype(int).values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
print(f'Train shape: {X_train.shape}')
print(f'Test shape: {X_test.shape}')
print('\nTraining Random Forest Classifier to distinguish condition (control vs activity) in OPCs...')
classifier = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
classifier.fit(X_train, y_train)
y_pred = classifier.predict(X_test)
print('\nClassification Report (on held-out test cells):')
print(classification_report(y_test, y_pred, target_names=['control', 'activity']))
print(f'Overall Accuracy: {accuracy_score(y_test, y_pred):.4f}')
print('\nCalculating SHAP values using TreeExplainer...')
explainer = shap.TreeExplainer(classifier)
shap_values_output = explainer.shap_values(X_test)
if isinstance(shap_values_output, list):
    shap_values = shap_values_output[1]
elif len(shap_values_output.shape) == 3:
    shap_values = shap_values_output[:, :, 1]
else:
    shap_values = shap_values_output
mean_abs_shap = np.abs(shap_values).mean(axis=0)
control_cells_mask = (adata_opc.obs['condition'] == 'control').values
activity_cells_mask = (adata_opc.obs['condition'] == 'activity').values
mean_expression_control = X[control_cells_mask].mean(axis=0)
mean_expression_activity = X[activity_cells_mask].mean(axis=0)
regulation_direction = np.where(mean_expression_activity > mean_expression_control, 'upregulated_in_activity', 'upregulated_in_control')
shap_df = pd.DataFrame({'gene_symbol': feature_gene_names, 'mean_abs_shap': mean_abs_shap, 'mean_expr_control': np.round(mean_expression_control, 4), 'mean_expr_activity': np.round(mean_expression_activity, 4), 'direction': regulation_direction}).sort_values('mean_abs_shap', ascending=False)
shap_df.to_csv(OUTPUT_CSV, index=False)
print(f'\nSaved all SHAP candidate genes to: {OUTPUT_CSV}')
print('\nTop 20 SHAP candidate genes in OPC-like cells with regulation details:')
print(shap_df.head(20).to_string(index=False))
print('\nCreating SHAP summary plot...')
plt.figure(figsize=(10, 6))
shap.summary_plot(shap_values, X_test, feature_names=feature_gene_names, max_display=20, show=False)
plt.title('Top 20 Genes Distinguishing Activity vs Control in OPC-like Cells (SHAP)', fontsize=14, pad=15)
plt.tight_layout()
plt.savefig(PLOT_FILE, dpi=300, bbox_inches='tight')
plt.close()
print(f'Saved SHAP summary plot to: {PLOT_FILE}')
print('\nDone.')
