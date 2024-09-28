import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
import matplotlib.pyplot as plt

# Function to compute molecular descriptors
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    descriptors = {
        "Molecular Weight": Descriptors.ExactMolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "Number of Hydrogen Donors": Descriptors.NumHDonors(mol),
        "Number of Hydrogen Acceptors": Descriptors.NumHAcceptors(mol),
        "Number of Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Aromatic Proportion": sum(atom.GetIsAromatic() for atom in mol.GetAtoms()) / mol.GetNumAtoms()
    }
    return descriptors

# Function to predict solubility (dummy function for demonstration)
def predict_solubility(smiles):
    return -1.0  # Dummy value, replace with actual prediction model

# Function to display molecular structure
def display_molecular_structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
        st.image(img, use_column_width=True)

# Title and instructions
st.title("Solubility Predictor")
st.write("Enter the SMILES notation of the molecule to predict its solubility and compute molecular descriptors.")

# Input SMILES string
smiles_input = st.text_area("SMILES Input (separate multiple SMILES with a newline):")
smiles_list = [smi.strip() for smi in smiles_input.split("\n") if smi.strip()]

if smiles_list:
    results = []
    for smiles in smiles_list:
        descriptors = compute_descriptors(smiles)
        if descriptors:
            solubility = predict_solubility(smiles)
            descriptors["Predicted Solubility (LogS)"] = solubility
            results.append(descriptors)
        else:
            st.error(f"Invalid SMILES string: {smiles}")

    if results:
        df = pd.DataFrame(results)
        st.write(df)

        # Download results as CSV
        csv = df.to_csv(index=False)
        st.download_button("Download Results", csv, "results.csv", "text/csv")

        # Visualize molecular structures and solubility
        st.header("Molecular Structures and Predicted Solubility")
        for smiles, row in zip(smiles_list, results):
            st.subheader(f"SMILES: {smiles}")
            display_molecular_structure(smiles)
            st.write(f"Predicted Solubility (LogS): {row['Predicted Solubility (LogS)']}")
            st.write("---")

        # Plot predicted solubility
        st.header("Predicted Solubility Plot")
        fig, ax = plt.subplots()
        ax.bar(smiles_list, df["Predicted Solubility (LogS)"])
        ax.set_xlabel("SMILES")
        ax.set_ylabel("Predicted Solubility (LogS)")
        plt.xticks(rotation=90)
        st.pyplot(fig)
    else:
        st.write("No valid SMILES strings provided.")
else:
    st.write("Please enter one or more SMILES strings to get started.")

# Footer
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
            <p><a href='#'>Privacy Policy</a></p>
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
