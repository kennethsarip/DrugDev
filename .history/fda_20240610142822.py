import mols2grid
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt, MolLogP, NumHDonors, NumHAcceptors

st.title("FDA-Approved Drugs")
st.write("""
This section will provide an interactive catalog of FDA-approved drugs. You can browse and filter drugs by various molecular properties.
""")

@st.cache_data
def download_dataset():
    df = pd.read_csv(
        "https://www.cureffi.org/wp-content/uploads/2013/10/drugs.txt", sep="\t"
    ).dropna()
    return df

def calc_mw(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return ExactMolWt(mol)

def calc_logp(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return MolLogP(mol)

def calc_NumHDonors(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHDonors(mol)

def calc_NumHAcceptors(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHAcceptors(mol)

def calc_num_rotatable_bonds(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)

def calc_aromatic_proportion(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    aromatic_atoms = sum(atom.GetIsAromatic() for atom in mol.GetAtoms())
    return aromatic_atoms / mol.GetNumAtoms()

df = download_dataset().copy()

if 'generic_name' in df.columns:
    df.rename(columns={'generic_name': 'Name', 'smiles': 'SMILES', 'cns_drug': 'Central Nervous System Drug'}, inplace=True)

df["Molecular Weight (g/mol)"] = df.apply(lambda x: f"{calc_mw(x['SMILES']):.1f} g/mol", axis=1)
df["LogP"] = df.apply(lambda x: calc_logp(x["SMILES"]), axis=1)
df["Number of Hydrogen Donors"] = df.apply(lambda x: calc_NumHDonors(x["SMILES"]), axis=1)
df["Number of Hydrogen Acceptors"] = df.apply(lambda x: calc_NumHAcceptors(x["SMILES"]), axis=1)
df["Number of Rotatable Bonds"] = df.apply(lambda x: f"{calc_num_rotatable_bonds(x['SMILES'])}", axis=1)
df["Aromatic Proportion"] = df.apply(lambda x: f"{calc_aromatic_proportion(x['SMILES']):.3f}", axis=1)

st.header('Set Parameters')
st.write('*Note: Display compounds having values less than the following thresholds*')

weight_cutoff = st.slider(
    label="Molecular Weight (g/mol)",
    min_value=0,
    max_value=1000,
    value=500,
    step=10,
)
logp_cutoff = st.slider(
    label="LogP",
    min_value=-10,
    max_value=10,
    value=5,
    step=1,
)
NumHDonors_cutoff = st.slider(
    label="Number of Hydrogen Donors",
    min_value=0,
    max_value=15,
    value=5,
    step=1,
)
NumHAcceptors_cutoff = st.slider(
    label="Number of Hydrogen Acceptors",
    min_value=0,
    max_value=20,
    value=10,
    step=1,
)

df_result = df[df["Molecular Weight (g/mol)"].apply(lambda x: float(x.split()[0])) < weight_cutoff]
df_result2 = df_result[df_result["LogP"] < logp_cutoff]
df_result3 = df_result2[df_result2["Number of Hydrogen Donors"] < NumHDonors_cutoff]
df_result4 = df_result3[df_result3["Number of Hydrogen Acceptors"] < NumHAcceptors_cutoff]

df_result4.reset_index(drop=True, inplace=True)

st.write(f"The filtered dataset has {df_result4.shape[0]} rows and {df_result4.shape[1]} columns.")
st.write(df_result4)

required_columns = ["SMILES", "Name", "Molecular Weight (g/mol)", "LogP", "Number of Hydrogen Donors", "Number of Hydrogen Acceptors", "Number of Rotatable Bonds", "Aromatic Proportion", "Central Nervous System Drug"]
if all(col in df_result4.columns for col in required_columns):
    raw_html = mols2grid.display(df_result4,
                                subset=["img", "Name", "Molecular Weight (g/mol)", "LogP", "Number of Hydrogen Donors", "Number of Hydrogen Acceptors", "Number of Rotatable Bonds", "Aromatic Proportion", "Central Nervous System Drug"],
                                mapping={"smiles": "SMILES", "generic_name": "Name"})._repr_html_()
    components.html(raw_html, width=900, height=1100, scrolling=False)
else:
    st.error(f"Missing one of the required columns: {required_columns}")

# Change the slider color dynamically based on the selection
NB = weight_cutoff

col = f''' <style> div.stSlider > div[data-baseweb = "slider"] > div > div {{
    background: linear-gradient(to right, rgb(29, 63, 114) 0%, 
                                rgb(29, 63, 114) {NB}%, 
                                rgba(151, 166, 195, 0.25) {NB}%, 
                                rgba(151, 166, 195, 0.25) 100%); }} </style>'''

st.markdown(col, unsafe_allow_html=True)

st.write("---")
st.markdown("""
    <div style='display: flex; justify-content: space-around; padding: 20px; background-color: white;'>
        <div style='background-color: white; padding: 10px;'>
            <p><a href='#'>FAQ</a></p>
        </div>
        <div style='background-color: white; padding: 10px;'>
            <p><a href='#'>User Guide</a></p>
        </div>
        <div style='background-color: white; padding: 10px;'>
            <p><a href='#'>Private Policy</a></p>
        </div>
        <div style='background-color: white; padding: 10px;'>
            <p><a href='#'>Terms of Service</a></p>
        </div>
    </div>
    <hr style='margin: 20px 0; padding: 0;'>
    <div style='text-align: center; padding: 10px;'>
        &copy; 2024 Drug Development App. All rights reserved. | Contact: info@drugdevapp.com
    </div>
    """, unsafe_allow_html=True)
