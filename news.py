import streamlit as st
import requests
from datetime import datetime

# Define your News API key
api_key = 'c9b2ff3e3bc84129a63613e46623d585'

# Define a function to fetch the latest news articles
def fetch_latest_news():
    query = ('neurological+diseases+AND+(drug+discovery+OR+drug+development)')
    url = ('https://newsapi.org/v2/everything?'
           f'q={query}&'
           'sortBy=publishedAt&'
           f'apiKey={api_key}')
    response = requests.get(url)
    articles = response.json().get('articles', [])
    return articles

# Fetch the latest news articles
news_articles = fetch_latest_news()

# Title and introductory text for the news section
st.title("Latest News in Drug Discovery and Development for Neurological Diseases")
st.write("Stay updated with the latest advancements and research in the field of bioinformatics, drug discovery, and development for neurological diseases.")

# Display each news article
for article in news_articles:
    st.subheader(article["title"])
    if article.get("author"):
        st.write(f"**Author:** {article['author']}")
    if article.get("publishedAt"):
        date_published = datetime.strptime(article["publishedAt"], '%Y-%m-%dT%H:%M:%SZ')
        st.write(f"**Published on:** {date_published.strftime('%B %d, %Y')}")
    if article.get("urlToImage"):
        st.image(article["urlToImage"], use_column_width=True)
    if article.get("description"):
        st.write(f"**Summary:** {article['description']}")
    st.markdown(f"[Read more]({article['url']})")

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
