"""
01_benchmark_construction.py

Builds the labeled benchmark from the Weng et al. (2024) deep mutational
scanning landscape. Partitions KRAS variants into interface (switch I,
residues 25-40; switch II, residues 60-76) and distal (all other positions),
and labels each as a binding disruptor, non-disruptor, or excluded, based on
measured binding free energy change to RAF1.

Classical hydrolysis-blocking hotspots (G12, G13, Q61) are excluded from the
positive class, since their pathogenicity arises from impaired GTP hydrolysis
rather than a change in per-encounter binding affinity.

Input:  kras_step1_mutation_table.csv (raw Weng et al. measurements, one row
        per verified single missense variant)
Output: labeled benchmark with a `category` column:
        interface_disruptor | interface_nondisruptor |
        distal_disruptor    | distal_nondisruptor
"""
import pandas as pd

SWITCH_I = range(25, 41)
SWITCH_II = range(60, 77)
HYDROLYSIS_HOTSPOTS = {12, 13, 61}

DISRUPTOR_THRESHOLD = 0.5    # kcal/mol, weaker binding
NONDISRUPTOR_THRESHOLD = 0.25  # kcal/mol, magnitude


def region(position: int) -> str:
    if position in SWITCH_I or position in SWITCH_II:
        return "interface"
    return "distal"


def label(row) -> str:
    if row["position"] in HYDROLYSIS_HOTSPOTS:
        return "excluded_hotspot"
    if not row.get("conf_RAF1_RBD", True):
        return "excluded_low_confidence"
    ddg = row["ddG_RAF1_RBD"]
    if ddg > DISRUPTOR_THRESHOLD:
        return f"{region(row['position'])}_disruptor"
    if abs(ddg) < NONDISRUPTOR_THRESHOLD:
        return f"{region(row['position'])}_nondisruptor"
    return "excluded_ambiguous"


def main():
    df = pd.read_csv("../data/kras_step1_mutation_table.csv")
    df["category"] = df.apply(label, axis=1)

    print("Category counts:")
    print(df["category"].value_counts())

    df.to_csv("../data/kras_benchmark_labeled.csv", index=False)
    print("\nSaved: kras_benchmark_labeled.csv")


if __name__ == "__main__":
    main()
