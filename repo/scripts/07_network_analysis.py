"""
07_network_analysis.py

Builds a KRAS residue contact graph from the predicted wild-type complex
structure (edge between two residues if their C-alpha atoms are within
8.5 A) and computes the shortest-path distance, in graph steps, from each
distal allosteric pocket residue to the RAF1 binding interface (formally
defined as switch I, residues 25-40, and switch II, residues 60-76,
matching the interface definition used throughout the paper).

Requires: gemmi, networkx
    pip install gemmi networkx

Input: the AlphaFold 3 wild-type complex structure (mmCIF).
"""
import argparse

import gemmi
import networkx as nx
import numpy as np

AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

SWITCH_I = set(range(25, 41))
SWITCH_II = set(range(60, 77))
INTERFACE_RESIDUES = SWITCH_I | SWITCH_II

# distal allosteric hotspot residues identified by the pocket-mapping analysis
POCKET_RESIDUES = [16, 17, 18, 19, 20, 21, 57, 58, 59, 146]

EDGE_DISTANCE_CUTOFF = 8.5  # Angstrom, CA-CA


def build_contact_graph(cif_path: str, kras_chain: str = "A") -> nx.Graph:
    structure = gemmi.read_structure(cif_path)
    model = structure[0]

    ca_coords = {}
    for chain in model:
        if chain.name != kras_chain:
            continue
        for residue in chain:
            if residue.name not in AMINO_ACIDS:
                continue
            for atom in residue:
                if atom.name == "CA":
                    ca_coords[residue.seqid.num] = np.array(
                        [atom.pos.x, atom.pos.y, atom.pos.z]
                    )

    graph = nx.Graph()
    residues = list(ca_coords.items())
    for i in range(len(residues)):
        for j in range(i + 1, len(residues)):
            r1, p1 = residues[i]
            r2, p2 = residues[j]
            if np.linalg.norm(p1 - p2) < EDGE_DISTANCE_CUTOFF:
                graph.add_edge(r1, r2)
    return graph


def shortest_paths_to_interface(graph: nx.Graph):
    paths = {}
    for pocket_res in POCKET_RESIDUES:
        if pocket_res not in graph:
            continue
        best_path = None
        for interface_res in INTERFACE_RESIDUES:
            if interface_res in graph and nx.has_path(graph, pocket_res, interface_res):
                path = nx.shortest_path(graph, pocket_res, interface_res)
                if best_path is None or len(path) < len(best_path):
                    best_path = path
        if best_path is not None:
            paths[pocket_res] = best_path
    return paths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wt_cif", required=True,
                     help="Path to the AlphaFold 3 wild-type complex mmCIF")
    args = ap.parse_args()

    graph = build_contact_graph(args.wt_cif)
    paths = shortest_paths_to_interface(graph)

    steps = {res: len(path) - 1 for res, path in paths.items()}
    print("Shortest path (graph steps) from each pocket residue to the RAF1 interface:")
    for res, n_steps in sorted(steps.items()):
        print(f"  residue {res}: {n_steps} step(s)  path = {paths[res]}")

    values = list(steps.values())
    print(f"\nmin = {min(values)}, max = {max(values)}, mean = {np.mean(values):.2f}")


if __name__ == "__main__":
    main()
