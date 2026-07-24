#Author: Adam Zelezny
#Date: 2026-07-15

import os
import pandas as pd

def parse_gsea_csv(csv_path: str) -> list[tuple[str, float]]:
    if not os.path.exists(csv_path):
        print(f'Warning: File not found at {csv_path}')
        return []
    df = pd.read_csv(csv_path)
    cleaned_terms = []
    for term in df['Term']:
        if ' (GO:' in term:
            cleaned_terms.append(term.split(' (GO:')[0])
        else:
            cleaned_terms.append(term)
    df['Cleaned_Term'] = cleaned_terms
    df['Abs_NES'] = df['NES'].abs()
    df_sorted = df.sort_values(by='Abs_NES', ascending=False).head(20)
    results = []
    for _, row in df_sorted.iterrows():
        results.append((row['Cleaned_Term'], row['NES']))
    return results

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    gsea_dir = os.path.join(script_dir, 'gsea_results')
    comparisons = [('Malignant', 'R88Q vs H1047R', 'tumor_H1047R_vs_R88Q_control'), ('Malignant', 'R88Q vs C420R', 'tumor_C420R_vs_R88Q_control'), ('Malignant', 'H1047R vs C420R', 'tumor_H1047R_vs_C420R_direct'), ('Neurons', 'R88Q vs H1047R', 'neuron_H1047R_vs_R88Q_control'), ('Neurons', 'R88Q vs C420R', 'neuron_C420R_vs_R88Q_control'), ('Neurons', 'H1047R vs C420R', 'neuron_H1047R_vs_C420R_direct')]
    latex_lines = []
    latex_lines.append('\\documentclass[12pt]{article}')
    latex_lines.append('\\usepackage[utf8]{inputenc}')
    latex_lines.append('\\usepackage{geometry}')
    latex_lines.append('\\geometry{letterpaper, margin=1in}')
    latex_lines.append('\\begin{document}')
    latex_lines.append('')
    current_section = None
    for pop, subsec, folder in comparisons:
        csv_path = os.path.join(gsea_dir, folder, 'gseapy.gene_set.prerank.report.csv')
        top_20 = parse_gsea_csv(csv_path)
        if pop != current_section:
            if current_section is not None:
                latex_lines.append('')
            latex_lines.append(f'\\section{{{pop}}}')
            current_section = pop
        latex_lines.append(f'\\subsection{{{subsec}}}')
        latex_lines.append('\\begin{itemize}')
        for term, nes in top_20:
            latex_lines.append(f'    \\item {term} {nes:.2f} NES')
        latex_lines.append('\\end{itemize}')
        latex_lines.append('')
    latex_lines.append('\\end{document}')
    out_path = os.path.join(script_dir, 'summary.txt')
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(latex_lines))
    print(f'LaTeX summary successfully written to {out_path}!')
if __name__ == '__main__':
    main()
