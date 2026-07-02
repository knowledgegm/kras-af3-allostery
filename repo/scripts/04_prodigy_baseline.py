"""
04_prodigy_baseline.py

Computes the PRODIGY (contact-based) structure-based binding affinity
baseline, for direct comparison against AlphaFold 3 interface confidence
on the same distal and interface variant sets.

Requires PRODIGY: pip install prodigy-prot

PRODIGY is run on the same AlphaFold-3-predicted mutant structures used
for the EvoEF2 baseline (03_evoef2_baseline.py), so both structure-based
baselines are evaluated on identical structural input.
"""
import argparse
import glob
import re
import subprocess

import pandas as pd


def prodigy_dg(pdb_path: str, chains=("A", "B")) -> float | None:
    result = subprocess.run(
        ["prodigy", pdb_path, "--selection", chains[0], chains[1], "-q"],
        capture_output=True, text=True,
    )
    match = re.search(r"(-?\d+\.\d+)", result.stdout)
    return float(match.group(1)) if match else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wt_pdb", default="wt_complex_Repair.pdb")
    ap.add_argument("--model_glob", default="wt_complex_Repair_Model_*.pdb",
                     help="Glob for the AF3-predicted mutant structures "
                          "built by the EvoEF2 BuildMutant step")
    ap.add_argument("--mutations_csv", default="../data/kras_step3_descriptors_5model.csv")
    ap.add_argument("--out", default="../data/kras_prodigy_baseline_reproduced.csv")
    args = ap.parse_args()

    muts = pd.read_csv(args.mutations_csv)["mutation"].tolist()
    wt_dg = prodigy_dg(args.wt_pdb)
    print(f"Wild-type PRODIGY dG: {wt_dg}")

    rows = []
    model_files = sorted(glob.glob(args.model_glob))
    for mutation, model_pdb in zip(muts, model_files):
        dg = prodigy_dg(model_pdb)
        if dg is not None:
            rows.append({"mutation": mutation, "prodigy_ddG": round(dg - wt_dg, 3)})

    pd.DataFrame(rows).to_csv(args.out, index=False)
    print(f"Saved {len(rows)} PRODIGY ddG values -> {args.out}")


if __name__ == "__main__":
    main()
