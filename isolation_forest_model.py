import pandas as pd
import pickle
import argparse

from table_preprocessing import preprocess_tables

from sklearn.ensemble import IsolationForest
from sklearn.metrics import pairwise_distances


def train_isolation_forest(
    X,
    contamination = 0.1,
    n_estimators = 100,
    random_state = 42,
    **kwargs,
):
    model = IsolationForest(
        n_estimators=n_estimators,
        contamination=contamination,
        random_state=random_state,
        **kwargs,
    )
    model.fit(X)
    return model


def predict_binary_vector(model, X):
    raw_predictions = model.predict(X)
    binary_vector = (raw_predictions == 1).astype(int)
    return binary_vector


def anomaly_scores(model, X):
    return model.score_samples(X)


def train_and_predict(
    X,
    contamination = 0.1,
    n_estimators = 100,
    random_state = 42,
):
    model = train_isolation_forest(
        X,
        contamination=contamination,
        n_estimators=n_estimators,
        random_state=random_state,
    )
    binary_vector = predict_binary_vector(model, X)
    scores = anomaly_scores(model, X)
    return model, binary_vector, scores


def save_model(model, path):
    """Serialize the fitted model to a pickle file."""
    with open(path, "wb") as f:
        pickle.dump(model, f)
    print(f"Model saved to '{path}'")


def load_model(path):
    """Load a previously serialized IsolationForest model."""
    with open(path, "rb") as f:
        model = pickle.load(f)
    return model


def build_results_dataframe(
    ids,
    binary_vector,
    scores,
):
    df = pd.DataFrame(
        {
            "id": ids,
            "binary": binary_vector,
            "anomaly_score": scores,
        }
    )
    df = df.sort_values("anomaly_score").reset_index(drop=True)
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Train IsolationForest on 3 VCF files and output "
            "a binary vector per variant (CHROM:POS)."
        )
    )
    parser.add_argument("table1", help="Path to first VCF file")
    parser.add_argument("table2", help="Path to second VCF file")
    parser.add_argument("table3", help="Path to third VCF file")
    parser.add_argument(
        "--contamination",
        type=float,
        default=0.1,
        help="Expected fraction of anomalous IDs (default: 0.1)",
    )
    parser.add_argument(
        "--n-estimators",
        type=int,
        default=100,
        help="Number of trees in the IsolationForest (default: 100)",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )
    parser.add_argument(
        "--no-scale",
        action="store_true",
        help="Skip StandardScaler during preprocessing",
    )
    parser.add_argument(
        "--output",
        default="binary_vector.csv",
        help="Output CSV file with id, binary label, and anomaly score (default: binary_vector.csv)",
    )
    parser.add_argument(
        "--save-model",
        default=None,
        help="If provided, save the fitted model to this .pkl path",
    )
    args = parser.parse_args()

    print("── Preprocessing tables ──────────────────────────────────────────")
    X, ids, feature_names, _ = preprocess_tables(
        args.table1,
        args.table2,
        args.table3,
        scale=not args.no_scale,
    )
    print(f"  Combined matrix shape : {X.shape}")
    print(f"  Number of IDs         : {len(ids)}")
    print(f"  Features per view     : {len(feature_names)}")

    print("\n── Training IsolationForest ──────────────────────────────────────")
    model, binary_vector, scores = train_and_predict(
        X,
        contamination=args.contamination,
        n_estimators=args.n_estimators,
        random_state=args.random_state,
    )
    n_anomalous = int((binary_vector == 0).sum())
    n_consistent = int((binary_vector == 1).sum())
    print(f"  Consistent IDs (1)  : {n_consistent}")
    print(f"  Anomalous  IDs (0)  : {n_anomalous}")

    print("\n── Saving results ────────────────────────────────────────────────")
    results_df = build_results_dataframe(ids, binary_vector, scores)
    results_df.to_csv(args.output, index=False)
    print(f"  Binary vector saved to '{args.output}'")

    if args.save_model:
        save_model(model, args.save_model)

    print("\nDone.")
