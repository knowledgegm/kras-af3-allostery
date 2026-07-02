#!/usr/bin/env python3
"""
run_admet.py  --  ADMET / developability panel for a list of molecules.
Runs on a normal laptop (CPU). Two layers:
  (1) RDKit rule-based panel (instant): physicochemical props, Lipinski/Veber/Ghose/Muegge,
      ESOL solubility, GI-absorption & BBB estimates, bioavailability score, SA, structural alerts.
  (2) OPTIONAL trained-model ADMET via admet-ai (toxicity, hERG, CYP, clearance, Caco-2, etc.)
      -> needs ~a few GB disk; install once with:  pip install admet-ai
Usage:
  pip install rdkit pandas
  python run_admet.py --input designs.csv --smiles_col smiles --out admet_out.csv
  (add --use_admet_ai to also run the heavy trained model if installed)
"""
import argparse, os, sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED, rdMolDescriptors as rdMD, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir,'SA_Score')); import sascorer

params=FilterCatalogParams()
for c in [FilterCatalogParams.FilterCatalogs.PAINS_A,FilterCatalogParams.FilterCatalogs.PAINS_B,
          FilterCatalogParams.FilterCatalogs.PAINS_C,FilterCatalogParams.FilterCatalogs.BRENK]:
    params.AddCatalog(c)
CAT=FilterCatalog.FilterCatalog(params)

def esol(m):
    mw=Descriptors.MolWt(m); logp=Crippen.MolLogP(m); rb=Descriptors.NumRotatableBonds(m)
    ap=sum(1 for a in m.GetAtoms() if a.GetIsAromatic())/max(m.GetNumHeavyAtoms(),1)
    return 0.16-0.63*logp-0.0062*mw+0.066*rb-0.74*ap

def rdkit_panel(smi):
    m=Chem.MolFromSmiles(smi)
    if m is None: return {}
    mw=Descriptors.MolWt(m); logp=Crippen.MolLogP(m); tpsa=rdMD.CalcTPSA(m)
    hbd=Lipinski.NumHDonors(m); hba=Lipinski.NumHAcceptors(m); rb=Descriptors.NumRotatableBonds(m)
    nar=rdMD.CalcNumAromaticRings(m); fsp3=rdMD.CalcFractionCSP3(m)
    lip=sum([mw>500,logp>5,hbd>5,hba>10])
    return dict(MW=round(mw),clogP=round(logp,2),TPSA=round(tpsa),HBD=hbd,HBA=hba,RotB=rb,
        ArRings=nar,Fsp3=round(fsp3,2),logS_ESOL=round(esol(m),2),
        GI_absorption='High' if (tpsa<=131.6 and -1<=logp<=5.88) else 'Low',
        BBB='Likely' if (tpsa<=79 and 0.4<=logp<=6) else 'Unlikely',
        Lipinski_viol=lip, Veber_pass=(rb<=10 and tpsa<=140),
        Bioavail_score=0.55 if (lip<=1 and tpsa<=140) else (0.17 if tpsa<=150 else 0.11),
        QED=round(QED.qed(m),2), SA=round(sascorer.calculateScore(m),1),
        struct_alerts=';'.join(sorted(set(e.GetDescription() for e in CAT.GetMatches(m))))[:60] or 'none')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True); ap.add_argument('--smiles_col',default='smiles')
    ap.add_argument('--out',default='admet_out.csv'); ap.add_argument('--use_admet_ai',action='store_true')
    a=ap.parse_args()
    df=pd.read_csv(a.input)
    rows=[{**{a.smiles_col:s}, **rdkit_panel(s)} for s in df[a.smiles_col]]
    out=pd.DataFrame(rows)
    if a.use_admet_ai:
        try:
            from admet_ai import ADMETModel
            mdl=ADMETModel()
            preds=mdl.predict(smiles=list(df[a.smiles_col]))
            out=pd.concat([out.reset_index(drop=True), preds.reset_index(drop=True)],axis=1)
            print("added trained-model ADMET endpoints")
        except Exception as e:
            print("admet-ai not available (pip install admet-ai). Skipping trained endpoints.", e)
    out.to_csv(a.out,index=False); print("wrote",a.out,"with",len(out),"molecules")

if __name__=='__main__': main()
