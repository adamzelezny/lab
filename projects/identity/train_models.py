import os
import sys
import numpy as np
from sklearn.metrics import accuracy_score, precision_recall_fscore_support

project_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(project_dir)

from src.models import (
    train_random_forest,
    train_mlp,
    evaluate_model,
    plot_confusion_matrix
)


def print_metric_explanations():
    print("=" * 75)
    print("METRIC EXPLANATIONS:")
    print("=" * 75)
    print("- Precision: Out of all cells PREDICTED as this type, what % were correct?")
    print("  (Higher precision = fewer false alarms)")
    print("- Recall: Out of all ACTUAL cells of this type, what % did the model find?")
    print("  (Higher recall = fewer missed cells)")
    print("- F1-Score: The balanced average of Precision and Recall.")
    print("  (The best single indicator of model quality, especially for rare cells)")
    print("- Support: The actual number of cells of this type in the test sample.")
    print("=" * 75)


def print_success_rates(y_true, y_pred):
    classes = sorted(list(np.unique(y_true)))
    precision, recall, f1, support = precision_recall_fscore_support(
        y_true, y_pred, labels=classes, zero_division=0
    )
    
    print(f"\n{'Cell Type':<25} | {'Precision':<10} | {'Recall':<10} | {'F1-Score':<10} | {'Support':<10}")
    print("-" * 75)
    for i, cls in enumerate(classes):
        print(f"{cls:<25} | {precision[i]:.4f}     | {recall[i]:.4f}  | {f1[i]:.4f}   | {support[i]:<10}")
    print("-" * 75)
    
    overall_acc = accuracy_score(y_true, y_pred)
    print(f"Overall Success Rate (Accuracy): {overall_acc:.4f}")
    print(f"Average F1-Score (Macro F1): {np.mean(f1):.4f}")


def main():
    processed_dir = os.path.join(project_dir, "data", "processed")
    train_path = os.path.join(processed_dir, "train_features.npz")
    test_path = os.path.join(processed_dir, "test_features.npz")
    
    output_dir = os.path.join(project_dir, "outputs")
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(train_path) or not os.path.exists(test_path):
        raise FileNotFoundError(
            f"Preprocessed features not found. Please run preprocess_dataset.py first."
        )
        
    train_data = np.load(train_path, allow_pickle=True)
    test_data = np.load(test_path, allow_pickle=True)
    
    X_train, y_train = train_data["X"], train_data["y"]
    X_test, y_test = test_data["X"], test_data["y"]
    
    class_names = sorted(list(np.unique(y_train)))
    
    print_metric_explanations()
    
    print("\nrunning random forest")
    rf_model = train_random_forest(
        X_train=X_train,
        y_train=y_train,
        random_state=42,
        n_estimators=100,
        class_weight=None
    )
    y_pred_rf, _ = evaluate_model(rf_model, X_test, y_test, class_names=class_names)
    print_success_rates(y_test, y_pred_rf)
    
    plot_path_rf = os.path.join(output_dir, "random_forest_standard_confusion_matrix.png")
    plot_confusion_matrix(
        y_true=y_test,
        y_pred=y_pred_rf,
        output_path=plot_path_rf,
        title="Random Forest Confusion Matrix (Standard)",
        class_names=class_names
    )
    
    print("\nrunning MLP")
    mlp_model = train_mlp(
        X_train=X_train,
        y_train=y_train,
        random_state=42,
        hidden_layer_sizes=(64, 32),
        alpha=0.0001,
        max_iter=1000
    )
    y_pred_mlp, _ = evaluate_model(mlp_model, X_test, y_test, class_names=class_names)
    print_success_rates(y_test, y_pred_mlp)
    
    plot_path_mlp = os.path.join(output_dir, "mlp_standard_confusion_matrix.png")
    plot_confusion_matrix(
        y_true=y_test,
        y_pred=y_pred_mlp,
        output_path=plot_path_mlp,
        title="MLP Confusion Matrix (Standard)",
        class_names=class_names
    )


if __name__ == "__main__":
    main()
