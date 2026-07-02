# AlphaFold 3 Interface Confidence Detects Distal Allosteric Modulation of KRAS-RAF1 Binding

Code and data for the paper *"AlphaFold 3 Interface Confidence Detects Distal
Allosteric Modulation of KRAS-RAF1 Binding That Structure-Based Energy
Functions Miss."*

This repository accompanies a submission to IEEE BIBM. It contains the
labeled benchmark, extracted AlphaFold 3 descriptors, structure-based
baseline results, and the analysis scripts needed to reproduce every
number reported in the paper.

## Summary

We test whether the interface confidence (ipTM) that AlphaFold 3 assigns to
a predicted KRAS-RAF1 complex can detect distal (non-interface) mutations
that allosterically weaken binding, using the deep mutational scanning
landscape of Weng et al. (*Nature*, 2024) as ground truth. The signal
recovers distal disruptors at an area under the curve of 0.79, is global
rather than local, and is not shared by two independent structure-based
energy functions (EvoEF2 and PRODIGY), which perform at chance on the same
task.

## Repository structure

```
.
├── data/
│   ├── kras_step1_mutation_table.csv        # full labeled benchmark, 3,453 variants
│   ├── kras_step3_descriptors_5model.csv    # AF3 ipTM/PAE descriptors, 77 complexes
│   ├── kras_step3_descriptors_all77.csv
│   ├── kras_baseline_features.csv           # simple physicochemical baseline features
│   ├── kras_evoef_baseline.csv              # EvoEF2 structure-based ddG results
│   ├── kras_prodigy_baseline.csv            # PRODIGY structure-based ddG results
│   ├── kras_ablation.csv                    # descriptor ablation results (Fig. 5)
│   ├── kras_alt_metrics.csv                 # alternative PAE-derived metrics
│   ├── kras_broken_contacts.csv             # per-residue contact-loss analysis (Fig. 7)
│   ├── kras_multieffector_generalization.csv
│   ├── kras_clinical_predictions_AF3.csv    # prospective clinical variant predictions
│   ├── kras_denovo_designs.csv              # de novo design candidates
│   ├── kras_designs_ADMET_trained.csv       # trained ADMET model results
│   ├── kras_top_designs_hardened.csv
│   ├── kras_network_steps.json              # pocket-to-interface path lengths
│   ├── kras_network_paths.json              # actual shortest-path residue chains
│   └── job_manifests/                       # AlphaFold 3 server job submission files
├── scripts/
│   ├── 01_benchmark_construction.py         # label variants: interface/distal x disruptor/non
│   ├── 02_extract_af3_descriptors.py        # parse AF3 outputs -> ipTM, PAE descriptors
│   ├── 03_evoef2_baseline.py                # EvoEF2 force-field ddG baseline
│   ├── 04_prodigy_baseline.py               # PRODIGY contact-based ddG baseline
│   ├── 05_statistical_analysis.py           # AUC, bootstrap CI, calibration, robustness
│   ├── 06_ablation_analysis.py              # which AF3 signal component matters (Fig. 5)
│   ├── 07_network_analysis.py               # pocket-to-interface contact-graph distance
│   ├── dock_into_kras_pocket.py             # AutoDock Vina docking pipeline
│   ├── run_admet.py                         # trained ADMET-AI profiling
│   └── run_kras_subsampling.py              # apo-ensemble pilot (early exploratory step)
└── requirements.txt
```

## Reproducing the core result

```bash
pip install -r requirements.txt

cd scripts
python 01_benchmark_construction.py
python 02_extract_af3_descriptors.py --af3_dir /path/to/af3_outputs
python 05_statistical_analysis.py
```

`05_statistical_analysis.py` reproduces the headline number (distal-task
AUC 0.79, 95% CI 0.67-0.89), the interface sanity check, the Spearman
correlation, the threshold-robustness sweep, the single-model-vs-averaged
comparison, and the calibration Brier score.

## Reproducing the structure-based baselines

The EvoEF2 and PRODIGY baselines require external tools:

```bash
# EvoEF2
git clone https://github.com/tommyhuangthu/EvoEF2.git
cd EvoEF2 && g++ -O3 -o EvoEF2 src/*.cpp

# PRODIGY
pip install prodigy-prot
```

Then, from a directory containing the repaired wild-type complex PDB:

```bash
python scripts/03_evoef2_baseline.py --evoef_bin ./EvoEF2/EvoEF2
python scripts/04_prodigy_baseline.py
```

## Reproducing the AlphaFold 3 predictions

AlphaFold 3 complex predictions were generated via the
[AlphaFold Server](https://alphafoldserver.com). The job submission
manifests in `data/job_manifests/` contain the exact sequences and ligand
specifications (KRAS variant, RAF1 RBD residues 56-131, GTP, Mg) used for
every prediction reported in the paper, and can be resubmitted directly.

## Notes on scope

This repository covers the core methodological pipeline: benchmark
construction, descriptor extraction, structure-based baselines, statistical
analysis, ablation, and network analysis. The exploratory translational
analyses described in the paper (druggable pocket mapping, drug repurposing
and de novo design, ADMET profiling, and prospective clinical variant
prediction) are included as data and, where a standalone script exists
(`dock_into_kras_pocket.py`, `run_admet.py`), as code; these are
explicitly labeled hypothesis-generating in the paper and are not required
to reproduce the central result.

## Citation

If you use this code or data, please cite:

```
[Author names], "AlphaFold 3 Interface Confidence Detects Distal Allosteric
Modulation of KRAS-RAF1 Binding That Structure-Based Energy Functions Miss,"
IEEE BIBM 2026 (submitted).
```

## License

[Choose a license before publishing, e.g., MIT for code, CC-BY-4.0 for data.]
