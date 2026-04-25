FROM python:3.10

WORKDIR /app

# Install regex directly
RUN pip install --no-cache-dir streamlit regex pandas biopython huggingface_hub python-dotenv scikit-learn numpy plotly bokeh dna_features_viewer matplotlib fastparquet pyarrow

# Copy everything
COPY . .

# Explicitly expose the port HF wants
EXPOSE 7860

# Run with the exact port HF expects for Docker
CMD ["streamlit", "run", "scripts/app.py", "--server.port=7860", "--server.address=0.0.0.0"]
