import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def load_data(file1, file2):
    matrix1 = np.loadtxt(file1, delimiter=",")
    matrix2 = np.loadtxt(file2, delimiter=",")
    return matrix1, matrix2

def perform_pca(matrix1, matrix2):
    X = np.vstack((matrix1, matrix2))
    labels = np.concatenate((np.zeros(matrix1.shape[0]), np.ones(matrix2.shape[0])))
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    return X_pca, labels, pca.explained_variance_ratio_

def plot_pca(X_pca, labels, explained_variance, save_path="prep/pca_plot.png"):
    plt.figure(figsize=(8, 6))
    
    plt.scatter(
        X_pca[labels == 0, 0],
        X_pca[labels == 0, 1],
        color="#3182bd",
        label="Pre-activation",
        alpha=0.8,
        edgecolors="k",
        s=80
    )
    
    plt.scatter(
        X_pca[labels == 1, 0],
        X_pca[labels == 1, 1],
        color="#e6550d",
        label="Activation",
        alpha=0.8,
        edgecolors="k",
        s=80
    )
    
    plt.xlabel(f"Principal Component 1 ({explained_variance[0]*100:.1f}% Variance)")
    plt.ylabel(f"Principal Component 2 ({explained_variance[1]*100:.1f}% Variance)")
    plt.title("PCA Dimensionality Reduction: Visualizing Dataset Differences")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    print(f"Successfully generated and saved PCA plot to {save_path}")

def main():
    matrix1, matrix2 = load_data("pre_activation.csv", "activation.csv")
    X_pca, labels, explained_variance = perform_pca(matrix1, matrix2)
    plot_pca(X_pca, labels, explained_variance)

if __name__ == "__main__":
    main()
