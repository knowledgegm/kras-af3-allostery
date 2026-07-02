"""
05_statistical_analysis.py

Runs the core statistical analysis of the paper:
  - AUC (distal and interface tasks) with bootstrap 95% confidence intervals
  - Spearman correlation between ipTM and measured binding change
  - Threshold-robustness sweep
  - Single-model vs. five-model-averaged comparison
  - Logistic calibration of ipTM to a probability of disruption,
    with cross-validated Brier score

Input: the labeled benchmark (01_benchmark_construction.py output) merged
with AF3 descriptors (02_extract_af3_descriptors.py output).
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss
from sklearn.model_selection import LeaveOneOut, cross_val_predict


def auc(scores: np.ndarray, labels: np.ndarray) -> float:
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    return float(np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg]))


def bootstrap_ci(scores: np.ndarray, labels: np.ndarray, n_boot: int = 5000, seed: int = 0):
    rng = np.random.default_rng(seed)
    idx = np.arange(len(scores))
    point = auc(scores, labels)
    boots = []
    for _ in range(n_boot):
        j = rng.choice(idx, len(idx), replace=True)
        if len(set(labels[j])) > 1:
            boots.append(auc(scores[j], labels[j]))
    return point, float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def main():
    df = pd.read_csv("../data/kras_step3_descriptors_5model.csv")

    # --- distal task (headline result) ---
    distal = df[df["category"].isin(["distal_disruptor", "distal_nondisruptor"])].copy()
    y = (distal["category"] == "distal_disruptor").astype(int).values
    scores = -distal["iptm5"].values  # lower ipTM = more disrupted

    point, lo, hi = bootstrap_ci(scores, y)
    print(f"Distal task AUC: {point:.3f}  95% CI [{lo:.3f}, {hi:.3f}]")

    rho, p = spearmanr(distal["iptm5"], distal["ddG"])
    print(f"Spearman rho (ipTM vs measured binding change): {rho:.3f} (p={p:.2e})")

    # --- interface sanity check ---
    interface = df[df["category"].isin(["interface_disruptor", "distal_nondisruptor"])].copy()
    yi = (interface["category"] == "interface_disruptor").astype(int).values
    print(f"Interface task AUC: {auc(-interface['iptm5'].values, yi):.3f}")

    # --- threshold robustness sweep ---
    print("\nThreshold robustness:")
    for cut in [0.4, 0.5, 0.6, 0.75, 1.0]:
        sub = df[(df["category"] == "distal_nondisruptor") |
                  ((df["category"] == "distal_disruptor") & (df["ddG"] > cut))]
        yy = (sub["category"] == "distal_disruptor").astype(int).values
        if len(set(yy)) > 1:
            print(f"  threshold {cut} kcal/mol -> AUC {auc(-sub['iptm5'].values, yy):.3f}")

    # --- single-model vs five-model averaging ---
    print(f"\nSingle-model ipTM AUC: {auc(-distal['iptm0'].values, y):.3f}")
    print(f"Five-model-averaged ipTM AUC: {auc(-distal['iptm5'].values, y):.3f}")

    # --- calibration ---
    X = distal[["iptm5"]].values
    proba = cross_val_predict(LogisticRegression(), X, y, cv=LeaveOneOut(), method="predict_proba")[:, 1]
    brier = brier_score_loss(y, proba)
    print(f"\nCross-validated Brier score: {brier:.3f}")


if __name__ == "__main__":
    main()
