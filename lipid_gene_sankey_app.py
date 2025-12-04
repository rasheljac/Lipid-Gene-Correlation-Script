import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from io import BytesIO, StringIO

st.set_page_config(page_title="Lipid-Gene Sankey Analysis", layout="wide")

st.title("Lipid-Gene Sankey Diagram Generator")
st.markdown("**White vs Beige Adipocyte Analysis**")

# Sidebar for file uploads
st.sidebar.header("Upload Data Files")

transcriptome_file = st.sidebar.file_uploader("Upload Transcriptome CSV", type=['csv'])
lipid_file = st.sidebar.file_uploader("Upload Lipid Data CSV", type=['csv'])

# Parameters
st.sidebar.header("Analysis Parameters")
gene_fc_threshold = st.sidebar.slider("Gene log2FC Threshold", 0.5, 3.0, 1.5, 0.1)
lipid_fc_threshold = st.sidebar.slider("Lipid log2FC Threshold", 0.1, 2.0, 0.8, 0.1)
top_genes_count = st.sidebar.slider("Top Genes to Display", 10, 100, 40, 5)

# Column name inputs
st.sidebar.header("Column Names")
beige_gene_prefix = st.sidebar.text_input("Beige Gene Column Prefix", "Hannah_Beige_")
white_gene_prefix = st.sidebar.text_input("White Gene Column Prefix", "Hannah_White_")
beige_lipid_prefix = st.sidebar.text_input("Beige Lipid Column Prefix", "Beige_")
white_lipid_prefix = st.sidebar.text_input("White Lipid Column Prefix", "White_")
gene_id_col = st.sidebar.text_input("Gene ID Column", "SampleID")
lipid_id_col = st.sidebar.text_input("Lipid ID Column", "Metabolite")

if transcriptome_file and lipid_file:
    # Load data
    transcriptome = pd.read_csv(transcriptome_file)
    lipid_data = pd.read_csv(lipid_file)
    
    # Clean transcriptome
    transcriptome = transcriptome[transcriptome[gene_id_col] != 'Class']
    transcriptome = transcriptome[transcriptome[gene_id_col].notna()]
    
    # Clean lipid data
    lipid_data = lipid_data[lipid_data[lipid_id_col] != 'Label']
    lipid_data = lipid_data[lipid_data[lipid_id_col] != 'Metabolite']
    lipid_data = lipid_data[lipid_data[lipid_id_col].notna()]
    
    # Find columns
    beige_gene_cols = [c for c in transcriptome.columns if c.startswith(beige_gene_prefix)]
    white_gene_cols = [c for c in transcriptome.columns if c.startswith(white_gene_prefix)]
    beige_lipid_cols = [c for c in lipid_data.columns if c.startswith(beige_lipid_prefix)]
    white_lipid_cols = [c for c in lipid_data.columns if c.startswith(white_lipid_prefix)]
    
    st.success("Files loaded successfully!")
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Genes", len(transcriptome))
        st.metric("Beige Gene Samples", len(beige_gene_cols))
        st.metric("White Gene Samples", len(white_gene_cols))
    with col2:
        st.metric("Lipids", len(lipid_data))
        st.metric("Beige Lipid Samples", len(beige_lipid_cols))
        st.metric("White Lipid Samples", len(white_lipid_cols))
    
    if st.button("Generate Sankey Diagram", type="primary"):
        with st.spinner("Computing fold changes..."):
            # Process genes
            for col in beige_gene_cols + white_gene_cols:
                transcriptome[col] = pd.to_numeric(transcriptome[col], errors='coerce')
            
            transcriptome['Beige_Mean'] = transcriptome[beige_gene_cols].mean(axis=1)
            transcriptome['White_Mean'] = transcriptome[white_gene_cols].mean(axis=1)
            transcriptome['Gene_log2FC'] = np.log2((transcriptome['White_Mean'] + 1) / (transcriptome['Beige_Mean'] + 1))
            
            # Process lipids
            for col in beige_lipid_cols + white_lipid_cols:
                lipid_data[col] = pd.to_numeric(lipid_data[col], errors='coerce')
            
            lipid_data['Beige_Mean'] = lipid_data[beige_lipid_cols].mean(axis=1)
            lipid_data['White_Mean'] = lipid_data[white_lipid_cols].mean(axis=1)
            lipid_data['Lipid_log2FC'] = np.log2((lipid_data['White_Mean'] + 1) / (lipid_data['Beige_Mean'] + 1))
            lipid_data['Lipid_Class'] = lipid_data[lipid_id_col].str.split(' ').str[0]
            
            # Filter significant changes
            sig_genes = transcriptome[abs(transcriptome['Gene_log2FC']) > gene_fc_threshold].copy()
            sig_lipids = lipid_data[abs(lipid_data['Lipid_log2FC']) > lipid_fc_threshold].copy()
            
            sig_genes['Direction'] = sig_genes['Gene_log2FC'].apply(lambda x: 'Upregulated_White' if x > 0 else 'Downregulated_White')
            sig_lipids['Direction'] = sig_lipids['Lipid_log2FC'].apply(lambda x: 'Upregulated_White' if x > 0 else 'Downregulated_White')
            
            # Select top genes
            genes_per_direction = top_genes_count // 2
            top_genes_up = sig_genes[sig_genes['Direction'] == 'Upregulated_White'].nlargest(genes_per_direction, 'Gene_log2FC')
            top_genes_down = sig_genes[sig_genes['Direction'] == 'Downregulated_White'].nsmallest(genes_per_direction, 'Gene_log2FC')
            top_genes = pd.concat([top_genes_up, top_genes_down])
            
            # Build flows
            flows = []
            for _, row in sig_lipids.iterrows():
                flows.append({'source': row['Lipid_Class'], 'target': row['Direction'], 'value': 1})
            
            for _, row in top_genes.iterrows():
                flows.append({'source': row['Direction'], 'target': row[gene_id_col], 'value': abs(row['Gene_log2FC'])})
            
            flows_df = pd.DataFrame(flows)
            flows_agg = flows_df.groupby(['source', 'target'])['value'].sum().reset_index()
            
            # Create nodes
            all_nodes = list(pd.concat([flows_agg['source'], flows_agg['target']]).unique())
            node_dict = {node: idx for idx, node in enumerate(all_nodes)}
            
            flows_agg['source_idx'] = flows_agg['source'].map(node_dict)
            flows_agg['target_idx'] = flows_agg['target'].map(node_dict)
            
            # Node colors
            node_colors = []
            for node in all_nodes:
                if 'Upregulated' in node:
                    node_colors.append('rgba(255, 99, 71, 0.8)')
                elif 'Downregulated' in node:
                    node_colors.append('rgba(70, 130, 180, 0.8)')
                elif node in ['PC', 'LPC', 'PE', 'LPE']:
                    node_colors.append('rgba(255, 165, 0, 0.8)')
                elif node in ['DG', 'TG']:
                    node_colors.append('rgba(34, 139, 34, 0.8)')
                elif node in ['SM', 'CAR']:
                    node_colors.append('rgba(148, 0, 211, 0.8)')
                else:
                    node_colors.append('rgba(128, 128, 128, 0.8)')
            
            # Create Sankey
            fig = go.Figure(data=[go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color='black', width=0.5),
                    label=all_nodes,
                    color=node_colors
                ),
                link=dict(
                    source=flows_agg['source_idx'],
                    target=flows_agg['target_idx'],
                    value=flows_agg['value']
                )
            )])
            
            fig.update_layout(
                title='Lipid-Gene Sankey: White vs Beige Adipocytes',
                font=dict(size=12),
                width=1400,
                height=900
            )
            
            st.plotly_chart(fig, width='stretch')
            
            # Summary statistics
            st.header("Analysis Summary")
            summary = {
                'Total_Genes': len(transcriptome),
                'Significant_Genes': len(sig_genes),
                'Genes_Up_White': len(sig_genes[sig_genes['Gene_log2FC'] > 0]),
                'Genes_Down_White': len(sig_genes[sig_genes['Gene_log2FC'] < 0]),
                'Total_Lipids': len(lipid_data),
                'Significant_Lipids': len(sig_lipids),
                'Lipids_Up_White': len(sig_lipids[sig_lipids['Lipid_log2FC'] > 0]),
                'Lipids_Down_White': len(sig_lipids[sig_lipids['Lipid_log2FC'] < 0])
            }
            summary_df = pd.DataFrame([summary]).T
            summary_df.columns = ['Count']
            st.dataframe(summary_df)
            
            # Export options
            st.header("Export Results")
            
            col1, col2 = st.columns(2)
            
            with col1:
                # HTML export - use to_html() method which returns string
                html_str = fig.to_html()
                st.download_button(
                    label="Download Sankey (HTML)",
                    data=html_str,
                    file_name="lipid_gene_sankey.html",
                    mime="text/html"
                )
            
            with col2:
                # CSV export
                csv_buffer = BytesIO()
                summary_df.to_csv(csv_buffer)
                csv_buffer.seek(0)
                st.download_button(
                    label="Download Summary (CSV)",
                    data=csv_buffer,
                    file_name="sankey_summary.csv",
                    mime="text/csv"
                )
            
            # Export flows
            flows_csv = BytesIO()
            flows_agg.to_csv(flows_csv, index=False)
            flows_csv.seek(0)
            st.download_button(
                label="Download Flow Data (CSV)",
                data=flows_csv,
                file_name="sankey_flows.csv",
                mime="text/csv"
            )

else:
    st.info("Upload both transcriptome and lipid data files to begin analysis")
    
    st.markdown("---")
    st.markdown("### Expected File Format")
    st.markdown("**Transcriptome CSV:**")
    st.markdown("- First column: Gene IDs")
    st.markdown("- Subsequent columns: Sample expression values (e.g., Hannah_Beige_1, Hannah_White_1)")
    st.markdown("")
    st.markdown("**Lipid CSV:**")
    st.markdown("- First column: Lipid names (e.g., PC 16:0)")
    st.markdown("- Subsequent columns: Sample intensity values (e.g., Beige_1, White_1)")

st.markdown("---")
st.markdown("**TeSlaa X Kapelczak Metabolomics**")
