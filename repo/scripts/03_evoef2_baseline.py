"""
03_evoef2_baseline.py

Computes the EvoEF2 (force-field) structure-based binding free energy
change baseline, for direct comparison against AlphaFold 3 interface
confidence on the same distal and interface variant sets.

Requires EvoEF2 compiled and available on PATH:
    git clone https://github.com/tommyhuangthu/EvoEF2.git
    cd EvoEF2 && g++ -O3 -o EvoEF2 src/*.cpp

Steps:
  1. Extract the wild-type KRAS-RAF1 complex from the AlphaFold 3
     wild-type prediction (chain A = KRAS, chain B = RAF1 RBD).
  2. Repair the structure with EvoEF2's RepairStructure command.
  3. Build each mutant with EvoEF2's BuildMutant command.
  4. Compute binding energy for wild type and every mutant with
     ComputeBinding, and take the difference (Eq. 3 in the paper):
         ddG = G_mut - G_wt
"""
import argparse
import re
import subprocess

import pandas as pd


def evoef_mutant_code(mutation: str) -> str:
    """Convert e.g. 'A146F' to EvoEF2's mutant-file format 'AA146F;'."""
    wt, pos, mt = mutation[0], mutation[1:-1], mutation[-1]
    return f"{wt}A{pos}{mt};"


def compute_binding(evoef_bin: str, pdb_path: str) -> float | None:
    result = subprocess.run(
        [evoef_bin, "--command=ComputeBinding", "--pdb", pdb_path],
        capture_output=True, text=True,
    )
    match = re.search(r"Total\s*=\s*(-?\d+\.\d+)", result.stdout)
    return float(match.group(1)) if match else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evoef_bin", default="./EvoEF2")
    ap.add_argument("--wt_pdb", default="wt_complex_Repair.pdb")
    ap.add_argument("--mutations_csv", default="../data/kras_step3_descriptors_5model.csv")
    ap.add_argument("--out", default="../data/kras_evoef2_baseline_reproduced.csv")
    args = ap.parse_args()

    muts = pd.read_csv(args.mutations_csv)["mutation"].tolist()

    mutant_file = "mutant_list.txt"
    with open(mutant_file, "w") as f:
        f.write("\n".join(evoef_mutant_code(m) for m in muts) + "\n")

    wt_energy = compute_binding(args.evoef_bin, args.wt_pdb)
    print(f"Wild-type binding energy: {wt_energy}")

    subprocess.run(
        [args.evoef_bin, "--command=BuildMutant", "--pdb", args.wt_pdb,
         "--mutant_file", mutant_file],
        check=True,
    )

    rows = []
    for i, mutation in enumerate(muts):
        model_pdb = f"wt_complex_Repair_Model_{i + 1:04d}.pdb"
        energy = compute_binding(args.evoef_bin, model_pdb)
        if energy is not None:
            rows.append({"mutation": mutation, "evoef_ddG": round(energy - wt_energy, 3)})

    pd.DataFrame(rows).to_csv(args.out, index=False)
    print(f"Saved {len(rows)} EvoEF2 ddG values -> {args.out}")


if __name__ == "__main__":
    main()
