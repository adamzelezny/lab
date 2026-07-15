import os
from typing import Tuple, List
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from xgboost import XGBClassifier
from sklearn.metrics import classification_report, confusion_matrix


def train_random_forest(
    X_train: np.ndarray,
    y_train: np.ndarray,
    random_state: int = 42,
    n_estimators: int = 100,
    class_weight: str = None
) -> RandomForestClassifier:
    model = RandomForestClassifier(
        n_estimators=n_estimators,
        class_weight=class_weight,
        random_state=random_state,
        n_jobs=-1
    )
    model.fit(X_train, y_train)
    return model


def train_logistic_regression(
    X_train: np.ndarray,
    y_train: np.ndarray,
    random_state: int = 42,
    C: float = 1.0,
    class_weight: str = "balanced",
    max_iter: int = 1000
) -> LogisticRegression:
    model = LogisticRegression(
        C=C,
        class_weight=class_weight,
        random_state=random_state,
        max_iter=max_iter
    )
    model.fit(X_train, y_train)
    return model


def train_mlp(
    X_train: np.ndarray,
    y_train: np.ndarray,
    random_state: int = 42,
    hidden_layer_sizes: Tuple = (64, 32),
    alpha: float = 0.0001,
    max_iter: int = 1000
) -> MLPClassifier:
    model = MLPClassifier(
        hidden_layer_sizes=hidden_layer_sizes,
        activation="relu",
        solver="adam",
        alpha=alpha,
        max_iter=max_iter,
        random_state=random_state,
        early_stopping=True,
        validation_fraction=0.1
    )
    model.fit(X_train, y_train)
    return model


def train_xgboost(
    X_train: np.ndarray,
    y_train: np.ndarray,
    random_state: int = 42,
    n_estimators: int = 100,
    max_depth: int = 6,
    learning_rate: float = 0.1
) -> XGBClassifier:
    model = XGBClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=learning_rate,
        random_state=random_state,
        n_jobs=-1,
        eval_metric="mlogloss"
    )
    X_train_clean = np.ascontiguousarray(X_train, dtype=np.float32)
    y_train_clean = np.ascontiguousarray(y_train, dtype=np.int32)
    model.fit(X_train_clean, y_train_clean)
    return model


def evaluate_model(
    model,
    X_test: np.ndarray,
    y_test: np.ndarray,
    class_names: List[str] = None
) -> Tuple[np.ndarray, np.ndarray]:
    X_test_clean = np.ascontiguousarray(X_test, dtype=np.float32)
    y_pred = model.predict(X_test_clean)
    labels = class_names if class_names is not None else np.unique(y_test)
    conf_mat = confusion_matrix(y_test, y_pred, labels=labels)
    return y_pred, conf_mat


def plot_confusion_matrix(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    output_path: str,
    title: str = "Confusion Matrix",
    class_names: List[str] = None
) -> None:
    labels = class_names if class_names is not None else np.unique(y_true)
    conf_mat = confusion_matrix(y_true, y_pred, labels=labels)
    
    row_sums = conf_mat.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    conf_mat_normalized = conf_mat.astype('float') / row_sums
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        conf_mat_normalized,
        annot=True,
        fmt=".2f",
        cmap="Blues",
        xticklabels=labels,
        yticklabels=labels,
        cbar=True,
        square=True
    )
    plt.title(title, fontsize=16, pad=15)
    plt.ylabel("True Label", fontsize=12)
    plt.xlabel("Predicted Label", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
