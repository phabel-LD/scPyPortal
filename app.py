from shiny import App, module, reactive, render, ui
#from shiny.express import input, session
import asyncio
import pandas as pd
import plotly.express as px
from pathlib import Path
import scanpy as sc

from config import Config
#from database import get_session
from utils.cache_manager import CacheManager
from utils.data_loader import DataLoader

# Import modules
from modules.preprocessing import create_preprocessing_ui, create_preprocessing_server
from modules.clustering import create_clustering_ui, create_clustering_server
from modules.differential_expression import create_de_ui, create_de_server
from modules.trajectory import create_trajectory_ui, create_trajectory_server
from modules.export import create_export_ui, create_export_server

# Initialize managers
cache_manager = CacheManager()
data_loader = DataLoader(cache_manager)

def create_sidebar():
    """Create the application sidebar"""
    return ui.sidebar(
        ui.h4("📊 scPyPortal"),
        ui.input_file(
            "file_upload", 
            "Upload Single-Cell Data",
            accept=[".h5ad", ".h5mu"],
            multiple=False,
            width="100%"
        ),
        ui.output_ui("dataset_info_panel"),
        ui.output_ui("qc_status"),
        ui.output_plot("qc_plots", height="450px"),
        ui.hr(),
        ui.output_ui("analysis_controls"),
        width=400,
        open="desktop"
    )

#@module.ui
def create_main_interface():
    """Create the main tabbed interface"""
    return ui.navset_pill(
        ui.nav_panel(
            "📈 Quality Control",
            create_preprocessing_ui("preprocessing")
        ),
        ui.nav_panel(
            "🔍 Clustering",
            create_clustering_ui("clustering")
        ),
        ui.nav_panel(
            "📊 Differential Expression",
            create_de_ui("differential_expression")
        ),
        ui.nav_panel(
            "🔄 Trajectory Analysis",
            create_trajectory_ui("trajectory")
        ),
        ui.nav_panel(
            "💾 Export Results",
            create_export_ui("export")
        ),
        id="main_tabs"
    )

app_ui = ui.page_sidebar(
    create_sidebar(),
    create_main_interface(),
    title="scPyPortal - Single-Cell Analysis Platform",
    fillable=True
)

def server(input, output, session):
    # Global reactive values
    current_adata = reactive.Value(None)
    current_dataset_hash = reactive.Value(None)
    is_processing = reactive.Value(False)
    qc_run = reactive.Value(False)

    # Initialize module servers
    preprocessing_module = create_preprocessing_server(
        "preprocessing", current_adata, is_processing
    )
    clustering_module = create_clustering_server(
        "clustering", current_adata, is_processing
    )
    de_module = create_de_server(
        "differential_expression", current_adata, is_processing
    )
    trajectory_module = create_trajectory_server(
        "trajectory", current_adata, is_processing
    )
    export_module = create_export_server(
        "export", current_adata, is_processing
    )
    
    @reactive.Effect
    @reactive.event(input.file_upload)
    async def handle_file_upload():
        """Handle file upload asynchronously"""
        is_processing.set(True)
        
        try:
            file_infos = input.file_upload()
            if not file_infos:
                return
            
            file_info = file_infos[0]
            file_path = Path(file_info["datapath"])
            
            # Cache and load dataset
            dataset_cache = await cache_manager.cache_dataset(
                file_path, file_info["name"]
            )
            
            # Load the cached dataset
            adata = data_loader.load_cached_dataset(dataset_cache.file_hash)
            current_adata.set(adata)
            current_dataset_hash.set(dataset_cache.file_hash)
            
            # Update session info
            session.send_custom_message("dataset_loaded", {
                "filename": file_info["name"],
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars
            })
            
        except Exception as e:
            session.notification_show(f"Error loading file: {str(e)}", type="error")
        finally:
            is_processing.set(False)
    
    @reactive.Effect
    @reactive.event(input.preprocess_data)
    def _handle_preprocess_button():
        #print("=== SIDEBAR PREPROCESS BUTTON CLICKED ===")
        #print(f"Current adata: {current_adata()}")
        
        if current_adata() is None:
            return
    
        is_processing.set(True)
        try:
            adata = current_adata().copy()
            print(f"Running QC on {adata.n_obs} cells, {adata.n_vars} genes")
        
            # Basic QC calculations
            adata.var['mt'] = adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        
            # Update the current_adata with QC results
            current_adata.set(adata)
            qc_run.set(True)
        
            print("QC completed successfully!")
        
        except Exception as e:
            print(f"QC error: {e}")
            import traceback
            traceback.print_exc()
        finally:
            is_processing.set(False)

    @output
    @render.ui
    def dataset_info_panel():
        """Display dataset information"""
        if current_adata() is None:
            return ui.p(
                "Please upload a .h5ad or .h5mu file to begin analysis.",
                class_="text-muted"
            )
        
        adata = current_adata()
        return ui.div(
            ui.h5("Dataset Info"),
            ui.p(f"📁 {adata.uns.get('original_filename', 'Unknown')}"),
            ui.p(f"🔬 Cells: {adata.n_obs:,}"),
            ui.p(f"🧬 Genes: {adata.n_vars:,}"),
            ui.hr(),
            ui.h6("Metadata Fields"),
            ui.p(f"Sample annotations: {len(adata.obs.columns)}"),
            ui.p(f"Gene annotations: {len(adata.var.columns)}"),
            class_="dataset-info"
        )
    
    @render.ui
    def qc_status():
        """Display QC status after processing"""
        adata = current_adata()
        if adata is None:
            return ui.p("")
    
        # Check if QC metrics exist
        if 'n_genes_by_counts' in adata.obs.columns:
            return ui.div(
                ui.h6("✅ QC Complete"),
                ui.p(f"Mean genes/cell: {adata.obs['n_genes_by_counts'].mean():.0f}"),
                ui.p(f"Mean counts/cell: {adata.obs['total_counts'].mean():.0f}"),
                style="color: green; padding: 10px; background: #e8f5e9; border-radius: 5px; margin-top: 10px;"
            )
        return ui.p("")
    
    @render.plot
    def qc_plots():
        """Display QC violin plots"""
        print("=== QC_PLOTS FUNCTION CALLED ===")
        print(f"qc_run value: {qc_run()}")
        
        if not qc_run():
            print("Returning None - qc not run")
            return None

        adata = current_adata()
        print(f"adata is None: {adata is None}")

        # Show plots if QC was run OR if data already has QC metrics
        if adata is None or 'n_genes_by_counts' not in adata.obs.columns:
            print("Returning None - no data or no qc columns")
            return None

        # Only require qc_run() if you want to force manual QC
        # For now, show plots if data has the metrics
        if not qc_run() and 'n_genes_by_counts' in adata.obs.columns:
            print("Data already has QC metrics, showing plots")

        print("Creating matplotlib plot...")
        import matplotlib.pyplot as plt
    
        fig, axes = plt.subplots(1, 3, figsize=(6, 4))
    
        # Genes per cell
        axes[0].hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black')
        axes[0].set_xlabel('Genes per Cell')
        axes[0].set_ylabel('Frequency')
        axes[0].set_title('Genes per Cell')
    
        # Total counts
        axes[1].hist(adata.obs['total_counts'], bins=50, edgecolor='black')
        axes[1].set_xlabel('Total Counts')
        axes[1].set_ylabel('Frequency')
        axes[1].set_title('UMI Counts per Cell')
    
        # Mitochondrial %
        if 'pct_counts_MT' in adata.obs.columns:
            axes[2].hist(adata.obs['pct_counts_MT'], bins=50, edgecolor='black')
            axes[2].set_xlabel('Mitochondrial %')
            axes[2].set_ylabel('Frequency')
            axes[2].set_title('Mitochondrial %')
    
        plt.tight_layout()
        print("Returning matplotlib figure")
        return fig

    @output
    @render.ui
    def analysis_controls():
        """Dynamic analysis controls"""
        if current_adata() is None:
            return ui.p("Upload data to enable analysis controls.")
        
        return ui.div(
            ui.h5("Analysis Controls"),
            ui.input_action_button(
                "preprocess_data",
                "Run Preprocessing",
                class_="btn-primary",
                width="100%"
            ),
            ui.input_action_button(
                "clear_cache",
                "Clear Cache",
                class_="btn-outline-secondary",
                width="100%"
            ),
            class_="analysis-controls"
        )
    
    @output
    @render.text
    def session_info():
        """Display session information"""
        if current_dataset_hash() is None:
            return "No active session"
        
        return f"Session: {current_dataset_hash()[:16]}... | Cells: {current_adata().n_obs}"

app = App(app_ui, server)