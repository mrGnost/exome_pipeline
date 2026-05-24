from isolation_forest_model import load_model, predict_binary_vector, build_results_dataframe, anomaly_scores
import argparse


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Predict on 3 VCF files and output "
            "a binary vector per variant (CHROM:POS)."
        )
    )
    parser.add_argument("table1", help="Path to first VCF file")
    parser.add_argument("table2", help="Path to second VCF file")
    parser.add_argument("table3", help="Path to third VCF file")
    parser.add_argument(
        "--model",
        default="model.pkl",
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
    args = parser.parse_args()

    from table_preprocessing import preprocess_tables

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

    model = load_model(args.model)
    binary_vector = predict_binary_vector(model, X)
    scores = anomaly_scores(model, X)
    n_anomalous = int((binary_vector == 0).sum())
    n_consistent = int((binary_vector == 1).sum())
    print(f"  Consistent IDs (1)  : {n_consistent}")
    print(f"  Anomalous  IDs (0)  : {n_anomalous}")

    print("\n── Saving results ────────────────────────────────────────────────")
    results_df = build_results_dataframe(ids, binary_vector, scores, feature_names)
    results_df.to_csv(args.output, index=False)
    print(f"  Binary vector saved to '{args.output}'")
    print("\nDone.")


if __name__ == "__main__":
    main()
