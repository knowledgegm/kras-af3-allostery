"""
02_extract_af3_descriptors.py

Extracts interface confidence descriptors from AlphaFold 3 server outputs
for each predicted KRAS-RAF1 complex (KRAS variant + RAF1 RBD residues
56-131 + GTP + Mg, five models per variant).

For each of the five models, reads:
  - fold_<name>_summary_confidences_<n>.json  -> ipTM, pTM
  - fold_<name>_full_data_<n>.json            -> full PAE matrix,
                                                  per-token chain IDs,
                                                  contact probabilities

Computes, per model and then averaged across the five models:
  - ipTM (primary descriptor used throughout the paper)
  - mean / min / 10th-percentile interface PAE
    (interface = cross-chain residue pairs with contact probability > 0.5)
  - TM-iface: a template-modeling-style interface score (Eq. 1 in the paper)
  - PAE-energy: a soft-minimum energy-style aggregation (Eq. 2 in the paper)
  - mutated-residue PAE: mean PAE between the mutated residue and all RAF1
    residues, used to test whether detection is local or global

Expects AF3 output folders named `cplx_<mutation>/` under --af3_dir, each
containing the summary_confidences and full_data JSON files described above.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

D0 = 10.0   # reference distance (A) for the TM-iface score, Eq. 1
TAU = 5.0   # temperature for the PAE-energy soft-minimum, Eq. 2


def tm_iface(pae_interface: np.ndarray) -> float:
    """Eq. 1: TM-iface = (1/N) * sum_i  1 / (1 + (PAE_i / d0)^2)"""
    return float(np.mean(1.0 / (1.0 + (pae_interface / D0) ** 2)))


def pae_energy(pae_interface: np.ndarray) -> float:
    """Eq. 2: PAE-energy = -tau * log( (1/N) * sum_i exp(-PAE_i / tau) )"""
    return float(-TAU * np.log(np.mean(np.exp(-pae_interface / TAU))))


def descriptors_for_model(summary_path: str, full_data_path: str, mutated_pos: int | None):
    summary = json.load(open(summary_path))
    full = json.load(open(full_data_path))

    pae = np.array(full["pae"])
    chain_ids = np.array(full["token_chain_ids"])
    contact_probs = np.array(full["contact_probs"])

    chains = list(dict.fromkeys(chain_ids))
    is_a = chain_ids == chains[0]
    is_b = chain_ids == chains[1]
    cross_mask = np.outer(is_a, is_b) | np.outer(is_b, is_a)
    contact_mask = cross_mask & (contact_probs > 0.5)
    interface_pae = pae[contact_mask] if contact_mask.sum() > 0 else pae[cross_mask]

    out = {
        "iptm": summary["iptm"],
        "ptm": summary.get("ptm"),
        "pae_mean": float(interface_pae.mean()),
        "pae_min": float(interface_pae.min()),
        "pae_p10": float(np.percentile(interface_pae, 10)),
        "tm_iface": tm_iface(interface_pae),
        "pae_energy": pae_energy(interface_pae),
    }

    if mutated_pos is not None:
        # mean PAE between the mutated residue (chain A) and all chain-B (RAF1) residues
        # token indices assumed 1:1 with residue numbering within chain A
        chain_a_idx = np.where(is_a)[0]
        if mutated_pos - 1 < len(chain_a_idx):
            row = chain_a_idx[mutated_pos - 1]
            out["mut_to_RAF_PAE"] = float(pae[row, is_b].mean())

    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--af3_dir", default="../data/af3_outputs",
                     help="Directory containing cplx_<mutation>/ subfolders")
    ap.add_argument("--out", default="../data/kras_af3_descriptors.csv")
    args = ap.parse_args()

    rows = []
    for folder in sorted(glob.glob(os.path.join(args.af3_dir, "cplx_*"))):
        mutation = os.path.basename(folder).replace("cplx_", "").upper()
        mutated_pos = None
        digits = "".join(c for c in mutation if c.isdigit())
        if digits:
            mutated_pos = int(digits)

        model_descs = []
        for n in range(5):
            summary_files = glob.glob(os.path.join(folder, f"*summary_confidences_{n}.json"))
            full_files = glob.glob(os.path.join(folder, f"*full_data_{n}.json"))
            if not summary_files or not full_files:
                continue
            model_descs.append(descriptors_for_model(summary_files[0], full_files[0], mutated_pos))

        if not model_descs:
            continue

        avg = {f"{k}5": float(np.mean([d[k] for d in model_descs if k in d]))
               for k in model_descs[0]}
        avg["mutation"] = mutation
        rows.append(avg)

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    print(f"Extracted descriptors for {len(df)} complexes -> {args.out}")


if __name__ == "__main__":
    main()
