import streamlit as st
from streamlit_option_menu import option_menu

# Define the color scheme
NAVY_COLOR = "#1D3F72"  # Navy blue from the FDA logo
WHITE_COLOR = "#FFFFFF"  # White

# Add custom CSS for the font and color scheme
st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Roboto:wght@400;700&display=swap');

    body {
        font-family: 'Helvetica', sans-serif;
    }

    div.stSlider > div[data-baseweb="slider"] > div[data-testid="stTickBar"] > div {
        background: rgb(1 1 1 / 0%);
    }

    div.stSlider > div[data-baseweb="slider"] > div > div > div[role="slider"]{
        background-color: #1D3F72; 
        box-shadow: rgb(14 38 74 / 20%) 0px 0px 0px 0.2rem;
    }

    div.stSlider > div[data-baseweb="slider"] > div > div > div > div {
        color: #1D3F72;
    }

    div.stSlider > div[data-baseweb="slider"] > div > div {
        background: linear-gradient(to right, rgb(1, 183, 158) 0%, 
                                    rgb(1, 183, 158) 50%, 
                                    rgba(151, 166, 195, 0.25) 50%, 
                                    rgba(151, 166, 195, 0.25) 100%);
    }
    </style>
    """, unsafe_allow_html=True)

st.image("drugdev.PNG", width = 150)

# Apply the color scheme to the option menu
page = option_menu(
    menu_title=None,
    options=["Home", "Predictor", "Drugs", "Diseases", "News"],
    icons=["house", "fa-flask", "capsule", "heart-pulse", "newspaper"],  # Appropriate icons
    menu_icon="cast",
    default_index=0,
    orientation='horizontal',
    styles={
        "container": {"padding": "0!important", "background-color": WHITE_COLOR},
        "icon": {"color": NAVY_COLOR, "font-size": "18px"},
        "nav-link": {
            "font-size": "16px", "text-align": "center", "margin": "0px", 
            "--hover-color": NAVY_COLOR, "color": NAVY_COLOR,
            "background-color": WHITE_COLOR, "font-family": "'Helvetica', sans-serif",
        },
        "nav-link-selected": {"background-color": NAVY_COLOR, "color": WHITE_COLOR},
        "icon-selected": {"color": WHITE_COLOR}  # Make the icon white when selected
    }
)

# Initialize session state for page if not already done
if 'page' not in st.session_state:
    st.session_state.page = "Home"

# Redirect to the appropriate page
if page == "Home":
    st.session_state.page = "Home"
    exec(open("home.py").read())
elif page == "Predictor":
    st.session_state.page = "Predictor"
    exec(open("predictor.py").read())
elif page == "Drugs":
    st.session_state.page = "Drugs"
    exec(open("fda.py").read())
elif page == "Diseases":
    st.session_state.page = "Diseases"
    exec(open("diseases.py").read())
elif page == "News":
    st.session_state.page = "News"
    exec(open("news.py").read())