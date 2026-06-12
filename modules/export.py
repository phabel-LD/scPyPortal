from shiny import module, reactive, render, ui
import io
from datetime import datetime
import scanpy as sc

@module.ui
def create_export_ui():
    return ui.div(
        ui.h3("Export Analysis Results"),
        ui.layout_sidebar(
            ui.sidebar(
                ui.h5("Export Options"),
                ui.input_text("export_filename", "Filename", 
                            value=f"scAnalysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}"),
                ui.p("This will export your processed data including:", style="margin-top: 20px;"),
                ui.tags.ul(
                    ui.tags.li("QC metrics"),
                    ui.tags.li("Clustering results (if computed)"),
                    ui.tags.li("UMAP/PCA embeddings"),
                    ui.tags.li("Differential expression results"),
                    ui.tags.li("Pseudotime (if computed)"),
                    ui.tags.li("All cell/gene annotations")
                ),
                ui.download_button("download_h5ad", "Download .h5ad File", class_="btn-primary"),
                width=350
            ),
            ui.div(
                ui.output_ui("export_info"),
                class_="main-content"
            )
        )
    )

@module.server
def create_export_server(input, output, session, current_adata, is_processing):
    
    @render.ui
    def export_info():
        adata = current_adata()
        
        if adata is None:
            return ui.div(
                ui.h4("No Data Loaded"),
                ui.p("Please upload and analyze data before exporting.")
            )
        
        # Summarize what will be exported
        info_items = []
        
        info_items.append(ui.p(f"📊 **Cells:** {adata.n_obs:,}"))
        info_items.append(ui.p(f"🧬 **Genes:** {adata.n_vars:,}"))
        
        # Check what analyses have been done
        analyses_done = []
        
        if 'n_genes_by_counts' in adata.obs.columns:
            analyses_done.append("✅ QC metrics")
        
        if 'X_pca' in adata.obsm:
            analyses_done.append(f"✅ PCA ({adata.obsm['X_pca'].shape[1]} components)")
        
        if 'X_umap' in adata.obsm:
            analyses_done.append("✅ UMAP")
        
        if 'leiden' in adata.obs.columns:
            n_clusters = len(adata.obs['leiden'].unique())
            analyses_done.append(f"✅ Leiden clustering ({n_clusters} clusters)")
        
        if 'dpt_pseudotime' in adata.obs.columns:
            analyses_done.append("✅ Pseudotime trajectory")
        
        if 'rank_genes_groups' in adata.uns:
            analyses_done.append("✅ Differential expression")
        
        # Display summary
        return ui.div(
            ui.h4("Export Summary"),
            ui.h5("Dataset Size:"),
            *info_items,
            ui.hr(),
            ui.h5("Analyses Included:"),
            ui.tags.ul(
                *[ui.tags.li(analysis) for analysis in analyses_done]
            ) if analyses_done else ui.p("No analyses computed yet", style="color: orange;"),
            ui.hr(),
            ui.h5("File Format:"),
            ui.p("HDF5-based AnnData format (.h5ad) - compatible with Scanpy, Seurat, and other tools"),
            ui.hr(),
            ui.p("Click the download button above to save your processed data.", 
                 style="font-style: italic; color: #666;")
        )
    
    @render.download(
        filename=lambda: f"{input.export_filename()}.h5ad"
    )
    def download_h5ad():
        print("\n" + "="*60)
        print("EXPORT DOWNLOAD REQUESTED")
        print("="*60)
        
        adata = current_adata()
        
        if adata is None:
            print("❌ No data to export")
            return
        
        print(f"📦 Exporting AnnData object:")
        print(f"   - Cells: {adata.n_obs:,}")
        print(f"   - Genes: {adata.n_vars:,}")
        print(f"   - Observations: {list(adata.obs.columns)}")
        print(f"   - Embeddings: {list(adata.obsm.keys())}")
        
         # Write to temporary file first, then read as bytes
        import tempfile
        import os

        with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp:
            tmp_path = tmp.name

        try:
            # Write to temp file
            adata.write_h5ad(tmp_path)

            # Read the file as bytes
            with open(tmp_path, 'rb') as f:
                data = f.read()
            
            print("✅ Export complete!")
            print("="*60 + "\n")

            yield data

        finally:
            # Clean up temp file
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
    
    return {}