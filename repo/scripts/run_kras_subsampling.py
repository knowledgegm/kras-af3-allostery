#!/usr/bin/env python3
"""
run_kras_subsampling.py
Generate conformational ensembles for KRAS wild type and mutants by
AlphaFold2 MSA subsampling (ColabFold). Lower MSA depth + many seeds makes
AF2 sample alternative switch-I / switch-II conformations instead of one state.

Run this inside a Colab session (or any machine with a GPU) AFTER installing
localColabFold / ColabFold. See kras_step2_README.md for setup.

Design: two phases so the MMseqs2 MSA server is called only ONCE per sequence.
  Phase A: build one MSA (.a3m) per sequence with --msa-only
  Phase B: predict many times from that MSA, sweeping MSA depth and seeds
"""
import os, shutil, glob, subprocess, csv, sys

# ----------------------------- CONFIG -----------------------------
FASTA          = "kras_step2_sequences.fasta"   # input sequences (WT + variants)
OUTROOT        = "kras_ensembles"               # all output lands here
MSA_DIR        = os.path.join(OUTROOT, "_msas") # precomputed a3m files

# Subsampling sweep. Each (max_seq, max_extra_seq) is one MSA depth.
# Smaller depth = more conformational diversity. (del Alamo 2022; Wayment-Steele 2024)
DEPTHS         = [(16, 32), (32, 64), (64, 128)]   # add (256,512) for a deeper, more native-biased point
NUM_SEEDS      = 5        # random seeds per depth (more = bigger ensemble)
NUM_MODELS     = 1        # AF2 model params per seed (1 is fine with many seeds; up to 5 for more diversity)
NUM_RECYCLE    = 3
USE_DROPOUT    = True     # dropout at inference adds sampling diversity (keep on for ensembles)

# PILOT vs FULL: keep PILOT True for the first run to validate the pipeline cheaply.
PILOT          = True
PILOT_VARIANTS = ["WT_KRAS_1-169", "S17N", "D57G", "T2D", "G77A"]  # WT + 2 disruptors + 2 non-disruptors
# ------------------------------------------------------------------


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ""
        else:
            seqs[name] += line
    return seqs


def run(cmd):
    print("  $", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main():
    if not shutil.which("colabfold_batch"):
        sys.exit("colabfold_batch not found. Install ColabFold first (see README).")

    seqs = read_fasta(FASTA)
    if PILOT:
        seqs = {k: v for k, v in seqs.items() if k in PILOT_VARIANTS}
        print(f"PILOT mode: {len(seqs)} sequences -> {PILOT_VARIANTS}")
    os.makedirs(MSA_DIR, exist_ok=True)

    manifest = []

    for name, seq in seqs.items():
        var_dir = os.path.join(OUTROOT, name)
        os.makedirs(var_dir, exist_ok=True)

        # ---- Phase A: one MSA per sequence ----
        a3m = os.path.join(MSA_DIR, f"{name}.a3m")
        if not os.path.exists(a3m):
            tmp_in = os.path.join(MSA_DIR, f"{name}.fasta")
            with open(tmp_in, "w") as fh:
                fh.write(f">{name}\n{seq}\n")
            run(["colabfold_batch", "--msa-only", tmp_in, MSA_DIR])
            # colabfold writes <name>.a3m into MSA_DIR
        if not os.path.exists(a3m):
            # some versions nest the a3m; find it
            hits = glob.glob(os.path.join(MSA_DIR, f"{name}*.a3m"))
            if hits:
                a3m = hits[0]
            else:
                print(f"  ! MSA not found for {name}, skipping"); continue

        # ---- Phase B: predict across depths and seeds, reusing the MSA ----
        for (mseq, mextra) in DEPTHS:
            out_d = os.path.join(var_dir, f"d{mseq}_{mextra}")
            os.makedirs(out_d, exist_ok=True)
            cmd = ["colabfold_batch", a3m, out_d,
                   "--max-seq", str(mseq), "--max-extra-seq", str(mextra),
                   "--num-seeds", str(NUM_SEEDS), "--num-models", str(NUM_MODELS),
                   "--num-recycle", str(NUM_RECYCLE)]
            if USE_DROPOUT:
                cmd.append("--use-dropout")
            run(cmd)
            for pdb in glob.glob(os.path.join(out_d, "*.pdb")):
                manifest.append({"variant": name, "max_seq": mseq, "max_extra_seq": mextra,
                                 "pdb": os.path.relpath(pdb, OUTROOT)})

    with open(os.path.join(OUTROOT, "ensemble_manifest.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["variant", "max_seq", "max_extra_seq", "pdb"])
        wr.writeheader(); wr.writerows(manifest)

    print(f"\nDone. {len(manifest)} structures across {len(seqs)} variants.")
    print(f"Manifest: {os.path.join(OUTROOT, 'ensemble_manifest.csv')}")
    print("Hand kras_ensembles/ to Step 4 (descriptor extraction).")


if __name__ == "__main__":
    main()
