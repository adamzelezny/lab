import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report

def run_random_forest(X_train, X_test, y_train, y_test):
    print("--- Random Forest Classifier ---")
    model = RandomForestClassifier(random_state=42)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    print(classification_report(y_test, y_pred))
    return model

def run_neural_network(X_train, X_test, y_train, y_test):
    print("--- Neural Network (MLP) Classifier ---")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    #1 hidden layer, 16 neurons
    model = MLPClassifier(
        hidden_layer_sizes=(16,),
        activation='relu',
        solver='adam',
        max_iter=1000,
        random_state=42,
        early_stopping=True,
        validation_fraction=0.1
    )
    model.fit(X_train_scaled, y_train)
    y_pred = model.predict(X_test_scaled)
    print(classification_report(y_test, y_pred))
    return model

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(script_dir, "..", "..", "synthetic_mouse_brain_cells.csv")
    
    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"Dataset not found at '{csv_path}'. Please run labeled_data.py first to generate it."
        )
        
    df = pd.read_csv(csv_path)

    X = df.drop(columns=["is_astrocyte"])
    y = df["is_astrocyte"]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    run_random_forest(X_train, X_test, y_train, y_test)

    run_neural_network(X_train, X_test, y_train, y_test)

if __name__ == "__main__":
    main()