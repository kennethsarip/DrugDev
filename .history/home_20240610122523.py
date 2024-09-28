import streamlit as st
import plotly.express as px
import pandas as pd

# Sample data for user activity
data = {
    'Date': pd.date_range(start='2024-06-01', periods=10, freq='D'),
    'Solubility Predictions': [5, 7, 6, 8, 7, 9, 6, 7, 8, 9],
    'FDA-Approved Drug Searches': [3, 2, 4, 3, 5, 4, 3, 4, 2, 3]
}

df = pd.DataFrame(data)

# Title and introductory text
st.title("Drug Development App")
st.write("This app aims to provide comprehensive tools and information to assist in drug development. Explore the latest news, predict molecular solubility, and interact with a catalog of FDA-approved drugs.")

# User activity section
st.header("User Activity")
df_melted = df.melt(id_vars=['Date'], value_vars=['Solubility Predictions', 'FDA-Approved Drug Searches'], 
                    var_name='Activity Type', value_name='Count')
fig = px.line(df_melted, x='Date', y='Count', color='Activity Type', title=None)
fig.update_layout(
    xaxis_title='Date',
    yaxis_title='Number of Activities',
    legend_title_text='Activity Type',
    hovermode='x unified'
)
st.plotly_chart(fig)

# Main Content Sections
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
