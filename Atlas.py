import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# ---------------------------------------------------------
# 1. Page Configuration & Custom CSS (The "Mesmerizing" UI)
# ---------------------------------------------------------
st.set_page_config(page_title="Toxo Phospho-Atlas", layout="wide", page_icon="🧬")

# Optional: Force a dark theme look using Plotly templates
CHART_THEME = "plotly_dark"
COLOR_SCALE = "magma" # Looks stunning for heatmaps

st.title("🧬 Toxoplasma Phosphoproteomic Atlas")
st.markdown("Explore high-confidence phosphorylation events, sequence motifs, and protein-specific sites.")

# ---------------------------------------------------------
# 2. File Upload & Data Processing
# ---------------------------------------------------------
with st.sidebar:
    st.header("📂 Data Input")
    uploaded_file = st.file_uploader("Upload Merged Data (.xlsx or .csv)", type=['xlsx', 'csv'])
    
    st.markdown("---")
    st.markdown("**UI Settings**")
    theme_choice = st.selectbox("Chart Theme", ["plotly_dark", "plotly_white", "seaborn"])
    CHART_THEME = theme_choice

if uploaded_file is not None:
    with st.spinner("Processing massive dataset..."):
        # Read file
        if uploaded_file.name.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
        else:
            df = pd.read_excel(uploaded_file)
            
        df['probability'] = pd.to_numeric(df['probability'], errors='coerce')
        
        # Filter for Gold class
        high_conf = df[df['probability'] > 0.75].copy()
        high_conf = high_conf.dropna(subset=['13_mer']) # Ensure peptides exist
        
        if high_conf.empty:
            st.error("No high-confidence sites found.")
            st.stop()

        # ---------------------------------------------------------
        # 3. Global Statistics
        # ---------------------------------------------------------
        total_proteins = high_conf['Protein_id'].nunique()
        total_sites = len(high_conf)
        avg_prob = high_conf['probability'].mean()
        
        col1, col2, col3 = st.columns(3)
        col1.metric("Unique Proteins Enriched", f"{total_proteins:,}")
        col2.metric("Total Gold Phosphosites", f"{total_sites:,}")
        col3.metric("Average Site Confidence", f"{avg_prob:.3f}")
        st.markdown("---")

        # ---------------------------------------------------------
        # 4. Tabbed Interface for Detailed Views
        # ---------------------------------------------------------
        tab1, tab2, tab3 = st.tabs(["📊 Global Landscape", "🔥 Sequence Motif Heatmap", "🔍 Protein Explorer"])

        # --- TAB 1: Global Landscape ---
        with tab1:
            st.subheader("Top Hyper-phosphorylated Proteins")
            # Calculate top proteins
            protein_counts = high_conf['Protein_id'].value_counts().reset_index()
            protein_counts.columns = ['Protein_id', 'Phosphosite_Count']
            top_50 = protein_counts.head(50)
            
            fig_bar = px.bar(
                top_50, x='Protein_id', y='Phosphosite_Count', 
                color='Phosphosite_Count', color_continuous_scale='viridis',
                template=CHART_THEME
            )
            fig_bar.update_layout(height=400, margin=dict(l=0, r=0, b=0, t=30))
            st.plotly_chart(fig_bar, use_container_width=True)

        # --- TAB 2: Sequence Motif Heatmap ---
        with tab2:
            st.subheader("Amino Acid Frequency Around Phosphorylation Site")
            st.markdown("This heatmap reveals the 'motif'—the amino acids most commonly found near the phosphorylation site (Position 0).")
            
            # Generate positional matrix
            amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
            # Positions: -6 to +6
            positions = [f"P{i}" for i in range(-6, 7)] 
            
            # Initialize empty matrix
            matrix = pd.DataFrame(0, index=amino_acids, columns=positions)
            
            # Count amino acids
            valid_peps = [p for p in high_conf['13_mer'] if len(str(p)) == 13]
            for p in valid_peps:
                for i, aa in enumerate(p):
                    if aa in matrix.index:
                        matrix.iloc[matrix.index.get_loc(aa), i] += 1
            
            # Convert to percentages
            matrix_pct = matrix.div(matrix.sum(axis=0), axis=1) * 100
            
            fig_heat = px.imshow(
                matrix_pct, 
                labels=dict(x="Sequence Position (0 = Phosphosite)", y="Amino Acid", color="% Frequency"),
                x=positions, y=amino_acids,
                color_continuous_scale=COLOR_SCALE,
                aspect="auto",
                template=CHART_THEME
            )
            # Highlight the center (position 0)
            fig_heat.add_vline(x=6, line_width=2, line_dash="dash", line_color="white")
            fig_heat.update_layout(height=600)
            st.plotly_chart(fig_heat, use_container_width=True)

        # --- TAB 3: Single Protein Explorer ---
        with tab3:
            st.subheader("Deep Dive: Single Protein Analysis")
            
            prot_list = sorted(high_conf['Protein_id'].unique())
            selected_prot = st.selectbox("Select or type a Protein ID to investigate:", prot_list)
            
            if selected_prot:
                prot_data = high_conf[high_conf['Protein_id'] == selected_prot]
                
                c1, c2 = st.columns([1, 2])
                with c1:
                    st.metric(f"Phosphosites on {selected_prot}", len(prot_data))
                    
                    # Probability Distribution for this protein
                    fig_dist = px.histogram(
                        prot_data, x="probability", nbins=10, 
                        title="Site Confidence Distribution",
                        color_discrete_sequence=['#00FFCC'],
                        template=CHART_THEME
                    )
                    fig_dist.update_layout(height=300)
                    st.plotly_chart(fig_dist, use_container_width=True)
                    
                with c2:
                    st.markdown("**Identified Phosphopeptides (13-mers)**")
                    st.dataframe(
                        prot_data[['13_mer', 'probability']].sort_values('probability', ascending=False),
                        use_container_width=True,
                        height=400
                    )
else:
    # Beautiful landing state before file is uploaded
    st.info("👈 Please upload your merged dataset in the sidebar to generate the Atlas.")
    st.markdown("""
    ### What this tool will do:
    1. **Filter** your data strictly for high-confidence (Class I) phosphosites.
    2. **Visualize** the global distribution of hyper-phosphorylated proteins.
    3. **Generate** a positional heatmap showing kinase motif preferences.
    4. **Provide** an interactive database to search and isolate specific proteins of interest.
    """)