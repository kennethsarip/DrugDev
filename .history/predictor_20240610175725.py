import streamlit as st
import pandas as pd
import pickle
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

def AromaticProportion(m):
    aromatic_atoms = [m.GetAtomWithIdx(i).GetIsAromatic() for i in range(m.GetNumAtoms())]
    AromaticAtom = sum(aromatic_atoms)
    HeavyAtom = Descriptors.HeavyAtomCount(m)
    if HeavyAtom == 0:
        return 0
    return AromaticAtom / HeavyAtom

def generate(smiles):
    moldata = [Chem.MolFromSmiles(elem) for elem in smiles if Chem.MolFromSmiles(elem) is not None]
    
    descriptors = []
    for mol in moldata:
        desc_MolLogP = round(Descriptors.MolLogP(mol), 2)
        desc_MolWt = round(Descriptors.MolWt(mol), 1)
        desc_NumRotatableBonds = round(Descriptors.NumRotatableBonds(mol), 0)
        desc_AromaticProportion = round(AromaticProportion(mol), 3)

        descriptors.append([desc_MolLogP, desc_MolWt, desc_NumRotatableBonds, desc_AromaticProportion])

    columnNames = ["MolLogP", "MolWt", "NumRotatableBonds", "AromaticProportion"]
    return pd.DataFrame(descriptors, columns=columnNames)

st.title("Solubility Predictor")
st.write("""
Predict the solubility of molecules using SMILES notation.
""")

# User input section
st.header('User Input')
SMILES_input = st.text_area("Enter SMILES strings (one per line):", value="")
SMILES = SMILES_input.strip().split('\n')
SMILES = [s for s in SMILES if s]

if SMILES:
    valid_smiles = [s for s in SMILES if Chem.MolFromSmiles(s)]
    invalid_smiles = [s for s in SMILES if not Chem.MolFromSmiles(s)]
    
    if invalid_smiles:
        st.error(f"Invalid SMILES strings: {', '.join(invalid_smiles)}")
    
    if valid_smiles:
        st.header('Computed Molecular Descriptors')
        X = generate(valid_smiles)

        for index, row in X.iterrows():
            col1, col2 = st.columns([1, 4])
            with col1:
                mol = Chem.MolFromSmiles(valid_smiles[index])
                img = Draw.MolToImage(mol)
                st.image(img, use_column_width=True)
            with col2:
                st.markdown(f"""
                **Molecule {index + 1}:** {valid_smiles[index]}
                - **Molecular LogP:** {row['MolLogP']}
                - **Molecular Weight:** {row['MolWt']} g/mol
                - **Number of Rotatable Bonds:** {int(row['NumRotatableBonds'])}
                - **Aromatic Proportion:** {row['AromaticProportion']}
                """)

        load_model = pickle.load(open('solubility_model.pkl', 'rb'))
        prediction = load_model.predict(X)

        st.header('Predicted Solubility (LogS) Values')
        results = pd.DataFrame({'SMILES': valid_smiles, 'Predicted LogS': prediction.round(2)})
        st.write(results)

        csv = results.to_csv(index=False)
        st.download_button("Download Results", csv, "results.csv", "text/csv")

else:
    st.write("No valid SMILES input provided.")

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