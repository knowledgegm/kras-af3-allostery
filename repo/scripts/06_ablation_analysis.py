"""
06_ablation_analysis.py

Ablates which component of the AlphaFold 3 confidence output carries the
distal detection signal: the mutated-residue PAE (local), interface PAE
aggregations (mean, min, 10th-percentile, TM-iface, PAE-energy), ipTM
(single-model and five-model-averaged), and a combined logistic model
over all interface descriptors (leave-one-out cross-validated).

Input: kras_step3_descriptors_5model.csv (must include a `mut_to_RAF_PAE`
column, produced by 02_extract_af3_descriptors.py) merged with the
alternative-metric columns (kras_alt_metrics.csv).
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def auc(scores: np.ndarray, labels: np.ndarray) -> float:
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    return float(np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg]))


def best_direction_auc(df: pd.DataFrame, col: str, y: np.ndarray) -> float:
    s = df[col].values
    return max(auc(s, y), auc(-s, y))


def main():
    descriptors = pd.read_csv("../data/kras_step3_descriptors_5model.csv")
    alt_metrics = pd.read_csv("../data/kras_alt_metrics.csv")
    merged = descriptors.merge(alt_metrics, on=["mutation", "category"])

    distal = merged[merged["category"].str.startswith("distal")].copy()
    y = (distal["category"] == "distal_disruptor").astype(int).values

    components = {
        "Mutated-residue PAE (local)": "mut_to_RAF_PAE",
        "Min interface PAE": "pae_min",
        "Interface PAE energy (Eq. 2)": "pae_energy",
        "Mean interface PAE": "pae_mean",
        "10th-percentile interface PAE": "pae_p10",
        "TM-iface (Eq. 1)": "tm_max",
        "ipTM, single model": "iptm0",
        "ipTM, 5-model average": "iptm5",
    }

    print("Component-wise AUC on the distal task:")
    for label, col in components.items():
        if col in distal.columns:
            print(f"  {label:35s} AUC = {best_direction_auc(distal, col, y):.3f}")

    # combined logistic model over all interface descriptors, LOO-CV
    feature_cols = ["iptm5", "pae_mean", "pae_min", "pae_p10", "pae_energy", "tm_energy"]
    feature_cols = [c for c in feature_cols if c in distal.columns]
    X = distal[feature_cols].fillna(0).values
    proba = cross_val_predict(
        make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)),
        X, y, cv=LeaveOneOut(), method="predict_proba",
    )[:, 1]
    print(f"  {'All components combined (LOO-CV)':35s} AUC = {auc(proba, y):.3f}")


if __name__ == "__main__":
    main()
