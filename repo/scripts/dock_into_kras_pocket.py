#!/usr/bin/env python3
"""
dock_into_kras_pocket.py
Screen a library of small molecules into the KRAS allosteric pocket identified
in this project (the nucleotide-rim / SII-adjacent region where distal
RAF1-binding disruptors concentrate).

Proof-of-concept docking. Vina scores are approximate and the pocket here is in
the GTP-active state, so treat results as hypothesis-generating, not validated
hits. For a real screen, include known KRAS binders (e.g. sotorasib, adagrasib)
as positive controls, and consider modeling the inactive-state SII pocket too.

Inputs:
  - receptor PDBQT (prepared once with meeko: mk_prepare_receptor.py --read_pdb kras.pdb -o kras -p)
  - a CSV with columns: name, smiles  (your library, e.g. an FDA-approved set)
Outputs:
  - ranked CSV of best predicted affinity per molecule

Install: pip install rdkit meeko vina gemmi
"""
import argparse, contextlib, io
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina

# Pocket box from this project (centroid of the disruptor hotspot residues)
BOX_CENTER = [0.1, 1.5, -2.0]
BOX_SIZE   = [22, 22, 22]

def ligand_pdbqt(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError("bad SMILES")
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, AllChem.ETKDGv3()) != 0:
        raise ValueError("embed failed")
    AllChem.MMFFOptimizeMolecule(m)
    setups = MoleculePreparation().prepare(m)
    s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    return s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", default="kras_receptor.pdbqt")
    ap.add_argument("--library", required=True, help="CSV with columns name,smiles")
    ap.add_argument("--out", default="docking_results.csv")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    args = ap.parse_args()

    lib = pd.read_csv(args.library)
    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(args.receptor)
    v.compute_vina_maps(center=BOX_CENTER, box_size=BOX_SIZE)

    rows = []
    for r in lib.itertuples(index=False):
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                v.set_ligand_from_string(ligand_pdbqt(r.smiles))
                v.dock(exhaustiveness=args.exhaustiveness, n_poses=5)
                e = v.energies(n_poses=1)[0][0]
            rows.append({"name": r.name, "best_affinity_kcal_mol": round(float(e), 2)})
        except Exception as ex:
            rows.append({"name": r.name, "best_affinity_kcal_mol": f"fail: {str(ex)[:40]}"})

    out = pd.DataFrame(rows).sort_values("best_affinity_kcal_mol")
    out.to_csv(args.out, index=False)
    print(out.to_string(index=False))

if __name__ == "__main__":
    main()
