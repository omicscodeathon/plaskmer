import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

def calculate_gc_content(sequence):
    """Calculates GC percentage from a nucleotide sequence."""
    if not sequence or pd.isna(sequence):
        return 0.0
    seq = sequence.upper()
    g = seq.count('G')
    c = seq.count('C')
    return ((g + c) / len(seq)) * 100 if len(seq) > 0 else 0.0

def prep_analytics_data(df):
    """Prepares the raw Parquet dataframe for fast analysis."""
    # Drop rows where Sequence is missing for accurate length/GC analysis
    df_clean = df.dropna(subset=['Sequence']).copy()
    
    # Calculate GC content on the fly!
    df_clean['GC_Content'] = df_clean['Sequence'].apply(calculate_gc_content)
    return df_clean

def plot_organism_heatmap(df, organism):
    """Generates a choropleth map of Africa showing sequencing density for an organism."""
    # Filter by organism and count records per country
    df_filtered = df[df['Organism'] == organism]
    country_counts = df_filtered['Country'].value_counts().reset_index()
    country_counts.columns = ['Country', 'Record_Count']
    
    fig = px.choropleth(
        country_counts,
        locations="Country",
        locationmode="country names",
        color="Record_Count",
        scope="africa",
        color_continuous_scale="Viridis",
        title=f"Geographical Distribution of {organism} Sequences in Africa",
        labels={'Record_Count': 'Number of Sequences'}
    )
    fig.update_layout(margin={"r":0,"t":40,"l":0,"b":0})
    return fig

def plot_type_distribution(df, organism):
    """Generates a donut chart showing what Types of sequences exist for an organism."""
    df_filtered = df[df['Organism'] == organism]
    type_counts = df_filtered['Type'].value_counts().reset_index()
    type_counts.columns = ['Sequence_Type', 'Count']
    
    fig = px.pie(
        type_counts, 
        values='Count', 
        names='Sequence_Type', 
        hole=0.4,
        title=f"Sequence Types for {organism}",
        color_discrete_sequence=px.colors.qualitative.Pastel
    )
    return fig

def plot_sequence_quality_metrics(df, organism):
    """Generates a dual boxplot for Sequence Lengths and GC Content."""
    df_filtered = df[df['Organism'] == organism]
    
    # Create a figure with secondary y-axis logic using subplots
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=1, cols=2, subplot_titles=("Sequence Length Distribution", "GC Content (%) Distribution"))
    
    fig.add_trace(go.Box(y=df_filtered['Length'], name="Length (bp)", marker_color='indianred'), row=1, col=1)
    fig.add_trace(go.Box(y=df_filtered['GC_Content'], name="GC Content (%)", marker_color='lightseagreen'), row=1, col=2)
    
    fig.update_layout(title_text=f"Quality & Genomic Metrics: {organism}", showlegend=False)
    return fig

def plot_cross_country_comparison(df, countries):
    """Generates a stacked bar chart comparing organisms across selected countries."""
    df_filtered = df[df['Country'].isin(countries)]
    grouped = df_filtered.groupby(['Country', 'Organism']).size().reset_index(name='Count')
    
    fig = px.bar(
        grouped, 
        x="Country", 
        y="Count", 
        color="Organism", 
        title="Cross-Country Organism Comparison",
        barmode="stack",
        text="Count"
    )
    return fig

def plot_specific_type_heatmap(df, sequence_type):
    """Generates a choropleth map specifically for Plasmids or mRNAs."""
    df_type = df[df['Type'].str.contains(sequence_type, case=False, na=False)]
    if df_type.empty:
        return None
        
    country_counts = df_type['Country'].value_counts().reset_index()
    country_counts.columns = ['Country', 'Record_Count']
    
    fig = px.choropleth(
        country_counts,
        locations="Country",
        locationmode="country names",
        color="Record_Count",
        scope="africa",
        color_continuous_scale="Plasma", # Different color scale to distinguish from the main map
        title=f"Geographical Distribution of {sequence_type}s in Africa"
    )
    fig.update_layout(margin={"r":0,"t":40,"l":0,"b":0})
    return fig

def plot_type_scatter_and_bar(df, sequence_type):
    """Generates a Host Organism bar chart and a Length vs GC scatter plot."""
    df_type = df[df['Type'].str.contains(sequence_type, case=False, na=False)]
    if df_type.empty:
        return None, None
        
    # 1. Bar Chart: Which organisms hold these sequences?
    org_counts = df_type['Organism'].value_counts().reset_index()
    org_counts.columns = ['Organism', 'Count']
    fig_bar = px.bar(
        org_counts, x="Organism", y="Count", color="Organism",
        title=f"Top Host Organisms for {sequence_type} Sequences"
    )
    
    # 2. Scatter Plot: Length vs GC Content
    fig_scatter = px.scatter(
        df_type, x="Length", y="GC_Content", color="Country", 
        hover_data=['Accession', 'Organism'],
        title=f"{sequence_type} Signature: Length vs GC Content (%)",
        labels={"Length": "Sequence Length (bp)", "GC_Content": "GC Content (%)"}
    )
    # Use log scale for length because plasmids/mRNAs can vary wildly in size
    fig_scatter.update_xaxes(type="log") 
    
    return fig_bar, fig_scatter
