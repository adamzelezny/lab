import csv
import numpy as np

def generate_baseline_expression(cells, genes, seed=42):
    np.random.seed(seed)    
    gene_means = np.random.uniform(1.0, 6.0, genes)
    gene_stds = np.random.uniform(0.2, 1.0, genes)
    
    pre = np.zeros((cells, genes))
    for g in range(genes):
        pre[:, g] = np.random.normal(
            gene_means[g],
            gene_stds[g],
            cells
        )
        
    pre = np.clip(pre, 0, None)
    return pre, gene_stds

def simulate_activation(pre, gene_stds, genes, cells, seed=42):
    activation = pre.copy()
    
    upregulated = np.random.choice(genes, 30, replace=False)
    
    remaining = list(set(range(genes)) - set(upregulated))
    downregulated = np.random.choice(remaining, 20, replace=False)
    
    remaining = list(set(remaining) - set(downregulated))
    heterogeneous = np.random.choice(remaining, 15, replace=False)
    
    for g in upregulated:
        activation[:, g] += np.random.normal(
            1.5,
            0.3,
            cells
        )
        
    for g in downregulated:
        activation[:, g] -= np.random.normal(
            1.2,
            0.2,
            cells
        )
        
    for g in heterogeneous:
        activation[:, g] += np.random.normal(
            0,
            gene_stds[g] * 1.5,
            cells
        )
        
    activation += np.random.normal(0, 0.15, activation.shape)
    activation = np.clip(activation, 0, None)
    
    return activation, upregulated, downregulated, heterogeneous

def save_datasets(pre, activation, file1="pre_activation.csv", file2="activation.csv"):
    with open(file1, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(np.round(pre, 4))
        
    with open(file2, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(np.round(activation, 4))
        
    print("Created:")
    print(f" - {file1}")
    print(f" - {file2}")

def main():
    cells = 200
    genes = 200
    seed = 42
    
    pre, gene_stds = generate_baseline_expression(cells, genes, seed)
    
    activation, upregulated, downregulated, heterogeneous = simulate_activation(
        pre, gene_stds, genes, cells, seed
    )
    
    save_datasets(pre, activation)
    
    print(f"Upregulated genes: {sorted(upregulated)}")
    print(f"Downregulated genes: {sorted(downregulated)}")
    print(f"Heterogeneous genes: {sorted(heterogeneous)}")

if __name__ == "__main__":
    main()