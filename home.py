import streamlit as st

st.write("---")

col1, col2, col3 = st.columns(3)

with col1:
    st.subheader("Latest News")
    st.write("Stay updated with the latest advancements and research in drug development. This section provides real-time news updates from reliable sources.")

with col2:
    st.subheader("Solubility Predictor")
    st.write("""
    The Solubility Predictor allows you to predict the solubility of molecules based on their SMILES notation.
    
    - **Input a SMILES string:** Enter the SMILES notation of the molecule you want to analyze.
    - **Get solubility predictions:** Receive an estimated solubility value along with a confidence score.
    - **Interactive results:** Visualize the molecular structure and understand the factors influencing solubility.
    """)

with col3:
    st.subheader("Interactive Catalog")
    st.write("""
    The Interactive Catalog features a comprehensive list of FDA-approved drugs, allowing you to:
    
    - **Browse molecular properties:** Explore detailed information on molecular weight, solubility (LogP), hydrogen bond donors, and acceptors.
    - **Filter molecules:** Use Lipinski's descriptors to filter and find molecules that meet specific criteria.
    - **Visualize structures:** View the chemical structure of each drug and understand its properties.
    """)

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
