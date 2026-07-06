import numpy as np
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from xgboost import XGBClassifier

def load_matrices(file1, file2):
    matrix1 = np.loadtxt(file1, delimiter=",")
    matrix2 = np.loadtxt(file2, delimiter=",")
    print(f"Matrix 1 Shape: {matrix1.shape}")
    print(f"Matrix 2 Shape: {matrix2.shape}\n")
    return matrix1, matrix2


#global stats

def global_statistics(matrix1, matrix2):
    print("=== Global Statistics ===")

    print(
        f"Matrix 1: Mean = {np.mean(matrix1):.4f}, "
        f"Median = {np.median(matrix1):.4f}, "
        f"Std = {np.std(matrix1):.4f}"
    )

    print(
        f"Matrix 2: Mean = {np.mean(matrix2):.4f}, "
        f"Median = {np.median(matrix2):.4f}, "
        f"Std = {np.std(matrix2):.4f}"
    )

    print()


#k-s test

def ks_test(matrix1, matrix2):

    print("=== Column-wise Statistical Tests (K-S Test) ===")

    significant_cols = []

    for col in range(matrix1.shape[1]):

        ks_stat, p_val = stats.ks_2samp(
            matrix1[:, col],
            matrix2[:, col]
        )

        if p_val < 0.05:
            significant_cols.append((col, p_val))

    if significant_cols:

        print(
            f"Found {len(significant_cols)} columns with "
            "significantly different distributions:"
        )

        for col, p in significant_cols:
            print(f"  - Column {col}: p-value = {p:.4f}")

    else:
        print("No significant columns found.")

    print()

    return significant_cols


#random forest

def random_forest(matrix1, matrix2):

    print("Random forest")

    X = np.vstack((matrix1, matrix2))

    y = np.concatenate((
        np.zeros(matrix1.shape[0]),
        np.ones(matrix2.shape[0])
    ))

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.5,
        random_state=42,
        stratify=y
    )

    clf = RandomForestClassifier(
        n_estimators=50,
        random_state=42
    )

    clf.fit(X_train, y_train)

    predictions = clf.predict(X_test)

    accuracy = accuracy_score(y_test, predictions)

    print(f"Accuracy: {accuracy*100:.1f}%")

    if accuracy <= 0.60:
        print("Conclusion: datasets are significantly similar.")
    elif accuracy <= 0.75:
        print("Conclusion: datasets are similar.")
    else:
        print("Conclusion: datasets are mostly different.")

    print()

    return clf


#xgboost

def xgboost(matrix1, matrix2):

    print("XGBoost")

    X = np.vstack((matrix1, matrix2))

    y = np.concatenate((
        np.zeros(matrix1.shape[0]),
        np.ones(matrix2.shape[0])
    ))

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.5,
        random_state=42,
        stratify=y
    )

    clf = XGBClassifier(
        n_estimators=50,
        random_state=42,
        eval_metric="logloss"
    )

    clf.fit(X_train, y_train)

    predictions = clf.predict(X_test)

    accuracy = accuracy_score(y_test, predictions)

    print(f"Accuracy: {accuracy*100:.1f}%")

    if accuracy <= 0.60:
        print("Conclusion: datasets are significantly similar.")
    elif accuracy <= 0.75:
        print("Conclusion: datasets are similar.")
    else:
        print("Conclusion: datasets are mostly different.")

    print()

    return clf


def main():

    matrix1, matrix2 = load_matrices(
        "pre_activation.csv",
        "activation.csv"
    )

    global_statistics(matrix1, matrix2)

    ks_test(matrix1, matrix2)

    random_forest(matrix1, matrix2)

    xgboost(matrix1, matrix2)


if __name__ == "__main__":
    main()