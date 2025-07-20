import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
import pandas as pd
import matplotlib.pyplot as plt

# ----- Branding -----
st.set_page_config(page_title="RADPT: Rajpal Agrawal Drug Prediction Tool", layout="wide")
st.title("🌟 RADPT: Rajpal Agrawal Drug Prediction Tool")
st.markdown("""
A free and open-source oral absorption and drug-likeness predictor for MD Pharmacology students and researchers.  
**Created by Dr. Sagar Agrawal in honor of Shri Rajpal Agrawal.**
""")

# --- Drug input section ---
st.header("Enter Drugs for Comparison")

num_drugs = st.number_input("How many drugs do you want to compare?", min_value=1, max_value=10, value=2)
drugs = {}
for i in range(num_drugs):
    cols = st.columns([1, 2])
    name = cols[0].text_input(f"Drug {i+1} Name", key=f"name{i}")
    smiles = cols[1].text_input(f"SMILES for {name or f'Drug {i+1}'}", key=f"smiles{i}")
    if name and smiles:
        drugs[name] = smiles

# --- Descriptor and rules function ---
def evaluate_rules(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rb = Lipinski.NumRotatableBonds(mol)
    lipinski = mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10
    pfizer = logp > 3 and tpsa < 75
    gsk = logp <= 2.5 and mw <= 400
    golden_triangle = 200 <= mw <= 500 and -0.4 <= logp <= 5.6
    veber = tpsa <= 140 and rb <= 10
    egan = logp <= 5.88 and tpsa <= 131.6
    bioavailability_score = 0.55 if lipinski and veber else 0.11
    return {
        "MW": round(mw, 2), "LogP": round(logp, 2), "HBD": hbd, "HBA": hba,
        "TPSA": round(tpsa, 2), "RB": rb,
        "Lipinski": "✅" if lipinski else "❌",
        "Pfizer": "✅" if pfizer else "❌",
        "GSK": "✅" if gsk else "❌",
        "Golden Triangle": "✅" if golden_triangle else "❌",
        "Veber": "✅" if veber else "❌",
        "Egan": "✅" if egan else "❌",
        "Bioavailability Score": bioavailability_score
    }

# --- Main Analysis ---
if drugs:
    results = []
    for name, smi in drugs.items():
        row = evaluate_rules(smi)
        if row:
            row["Drug"] = name
            row["SMILES"] = smi
            results.append(row)
    if results:
        df = pd.DataFrame(results)
        st.subheader("Comparison Table")
        st.dataframe(df, hide_index=True, use_container_width=True)

        # --- Bar chart visualization ---
        st.subheader("Drug Property Bar Charts (Green = Acceptable Range)")
        acceptable_ranges = {
            'MW': (0, 500),
            'LogP': (-2, 5),
            'TPSA': (0, 140),
            'HBD': (0, 5),
            'HBA': (0, 10),
            'RB': (0, 10)
        }
        for prop, (min_val, max_val) in acceptable_ranges.items():
            fig, ax = plt.subplots(figsize=(6, 3.5))
            ax.bar(df['Drug'], df[prop], color='skyblue', edgecolor='black')
            ax.axhspan(min_val, max_val, color='lightgreen', alpha=0.3, label='Acceptable Range')
            ax.set_title(f'{prop} Comparison')
            ax.set_ylabel(prop)
            ax.set_xlabel('Drug')
            ax.grid(axis='y', linestyle='--', alpha=0.7)
            ax.legend()
            st.pyplot(fig)

        st.info("For correct SMILES, use [PubChem](https://pubchem.ncbi.nlm.nih.gov/), ChEMBL, or similar chemical databases.")

# --- Footer and credits ---
st.markdown("---")
st.markdown(
    "© 2024 Dr. Sagar Agrawal. Open-source and free to use for the scientific community. "
    "For queries, contact [Sagar Agrawal](mailto:sagaragrawal1506@gmail.com)."
)
