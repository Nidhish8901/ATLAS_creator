import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# ---------------------------------------------------------
# 1. Page Configuration & Custom CSS
# ---------------------------------------------------------
st.set_page_config(page_title="Toxo Phospho-Atlas", layout="wide", page_icon="🧬")

st.title("🧬 Toxoplasma Phosphoproteomic Atlas")
st.markdown("Explore comprehensive phosphorylation events across Gold, Silver, and Bronze confidence thresholds.")

# ---------------------------------------------------------
# 2. File Upload & Data Processing
# ---------------------------------------------------------
with st.sidebar:
    st.header("📂 Data Input")
    uploaded_file = st.file_uploader("Upload Merged Data (.xlsx or .csv)", type=['xlsx', 'csv'])
    
    st.markdown("---")
    st.markdown("**UI Settings**")
    CHART_THEME = st.selectbox("Chart Theme", ["plotly_dark", "plotly_white", "seaborn"])
    COLOR_SCALE = st.selectbox("Heatmap Color", ["magma", "viridis", "plasma", "inferno"])

if uploaded_file is not None:
    with st.spinner("Classifying and processing dataset..."):
        # Read file
        if uploaded_file.name.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
        else:
            df = pd.read_excel(uploaded_file)
            
        df['probability'] = pd.to_numeric(df['probability'], errors='coerce')
        df = df.dropna(subset=['probability', '13_mer']) # Ensure valid data
        
        # --- Classification Engine ---
        def assign_class(prob):
            if prob > 0.75: return 'Gold (>0.75)'
            elif prob >= 0.50: return 'Silver (0.50 - 0.75)'
            else: return 'Bronze (<0.50)'
            
        df['Confidence_Class'] = df['probability'].apply(assign_class)
        
        # Color mapping for consistency across all charts
        class_colors = {
            'Gold (>0.75)': '#FFD700',   # Gold
            'Silver (0.50 - 0.75)': '#C0C0C0', # Silver
            'Bronze (<0.50)': '#CD7F32'  # Bronze
        }

        # ---------------------------------------------------------
        # 3. Global Statistics (Including Silver & Bronze)
        # ---------------------------------------------------------
        total_proteins = df['Protein_id'].nunique()
        gold_count = len(df[df['Confidence_Class'] == 'Gold (>0.75)'])
        silver_count = len(df[df['Confidence_Class'] == 'Silver (0.50 - 0.75)'])
        bronze_count = len(df[df['Confidence_Class'] == 'Bronze (<0.50)'])
        
        st.markdown("### 📈 Dataset Overview")
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Unique Proteins", f"{total_proteins:,}")
        col2.metric("🏆 Gold Sites", f"{gold_count:,}")
        col3.metric("🥈 Silver Sites", f"{silver_count:,}")
        col4.metric("🥉 Bronze Sites", f"{bronze_count:,}")
        st.markdown("---")

        # ---------------------------------------------------------
        # 4. Tabbed Interface for Detailed Views
        # ---------------------------------------------------------
        tab1, tab2, tab3 = st.tabs(["📊 Global Landscape", "🔥 Motif Heatmaps", "🔍 Protein Explorer"])

        # --- TAB 1: Global Landscape (Stacked Bar) ---
        with tab1:
            st.subheader("Top 50 Hyper-phosphorylated Proteins (By Confidence)")
            
            # Count sites per protein AND class
            protein_class_counts = df.groupby(['Protein_id', 'Confidence_Class']).size().reset_index(name='Count')
            
            # Get the top 50 proteins based on TOTAL sites
            top_proteins = df['Protein_id'].value_counts().head(50).index
            filtered_counts = protein_class_counts[protein_class_counts['Protein_id'].isin(top_proteins)]
            
            # Create a stacked bar chart
            fig_bar = px.bar(
                filtered_counts, x='Protein_id', y='Count', color='Confidence_Class',
                color_discrete_map=class_colors,
                template=CHART_THEME, barmode='stack',
                category_orders={"Confidence_Class": ["Gold (>0.75)", "Silver (0.50 - 0.75)", "Bronze (<0.50)"]}
            )
            fig_bar.update_layout(height=500, xaxis={'categoryorder':'total descending'}, margin=dict(t=30))
            st.plotly_chart(fig_bar, use_container_width=True)

        # --- TAB 2: Sequence Motif Heatmaps (Interactive Toggle) ---
        with tab2:
            st.subheader("Positional Amino Acid Frequency Heatmaps")
            
            # Let the user choose which heatmap to view
            selected_class = st.radio(
                "Select Confidence Level to Generate Heatmap:",
                ["Gold (>0.75)", "Silver (0.50 - 0.75)", "Bronze (<0.50)"],
                horizontal=True
            )
            
            # Add a biological disclaimer for Silver/Bronze
            if selected_class != "Gold (>0.75)":
                st.caption(f"⚠️ *Note: {selected_class.split(' ')[0]} sites have lower localization probability. The motif pattern here may be 'noisier' because the exact center (Position 0) might be shifted to a neighboring Serine/Threonine.*")
            
            # Filter data for the selected heatmap
            heatmap_data = df[df['Confidence_Class'] == selected_class]
            
            if heatmap_data.empty:
                st.warning(f"No {selected_class} sites found to generate a heatmap.")
            else:
                amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
                positions = [f"P{i}" for i in range(-6, 7)] 
                matrix = pd.DataFrame(0, index=amino_acids, columns=positions)
                
                valid_peps = [p for p in heatmap_data['13_mer'] if len(str(p)) == 13]
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
                    template=CHART_THEME,
                    title=f"Motif Heatmap: {selected_class} Sites (n={len(valid_peps)})"
                )
                fig_heat.add_vline(x=6, line_width=2, line_dash="dash", line_color="white")
                fig_heat.update_layout(height=600)
                st.plotly_chart(fig_heat, use_container_width=True)

        # --- TAB 3: Single Protein Explorer (All Classes) ---
        with tab3:
            st.subheader("Deep Dive: Single Protein Analysis")
            
            prot_list = sorted(df['Protein_id'].unique())
            selected_prot = st.selectbox("Search for a Protein ID:", prot_list)
            
            if selected_prot:
                prot_data = df[df['Protein_id'] == selected_prot]
                
                c1, c2 = st.columns([1, 2])
                with c1:
                    st.metric(f"Total Sites on {selected_prot}", len(prot_data))
                    
                    # Probability Distribution colored by class
                    fig_dist = px.histogram(
                        prot_data, x="probability", nbins=15, 
                        title="Confidence Distribution",
                        color="Confidence_Class",
                        color_discrete_map=class_colors,
                        template=CHART_THEME,
                        category_orders={"Confidence_Class": ["Gold (>0.75)", "Silver (0.50 - 0.75)", "Bronze (<0.50)"]}
                    )
                    fig_dist.update_layout(height=350, showlegend=False)
                    st.plotly_chart(fig_dist, use_container_width=True)
                    
                with c2:
                    st.markdown(f"**Identified Phosphopeptides for {selected_prot}**")
                    
                    # Style the dataframe based on confidence class
                    display_df = prot_data[['13_mer', 'probability', 'Confidence_Class']].sort_values('probability', ascending=False)
                    st.dataframe(display_df, use_container_width=True, height=350)
else:
    st.info("👈 Please upload your dataset in the sidebar to generate the Atlas.")
