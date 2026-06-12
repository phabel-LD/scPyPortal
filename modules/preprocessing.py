from shiny import module, reactive, render, ui
import scanpy as sc
import plotly.express as px

##@module.ui
#def create_preprocessing_ui(module_id: str):
@module.ui
def create_preprocessing_ui():
    return ui.div(
        ui.h3("Quality Control & Preprocessing"),
        ui.layout_sidebar(
            ui.sidebar(
                ui.h5("QC Parameters"),
                ui.input_numeric("min_genes", "Min Genes per Cell", value=200),
                ui.input_numeric("max_genes", "Max Genes per Cell", value=5000),
                ui.input_numeric("max_mito", "Max Mitochondrial %", value=10),
                ui.input_action_button("run_qc", "Run QC", class_="btn-primary"),
                width=300
            ),
            ui.div(
                ui.output_ui("qc_results_ui"),
                ui.output_plot("qc_plots_module", height="600px"),
                class_="main-content"
            )
        )
    )

@module.server
def create_preprocessing_server(input, output, session, current_adata, is_processing):
    qc_results = reactive.Value(None)
    
    @reactive.Effect
    @reactive.event(input.run_qc)
    def _run_qc():
        print("\n" + "="*60)
        print("MODULE QC BUTTON CLICKED")
        print("="*60)
        print(f"current_adata is None: {current_adata() is None}")  # Debug
    
        if current_adata() is None:
            print("-> No data loaded, returning")  # Debug
            return
    
        is_processing.set(True)
        print("-> Starting QC processing...")  # Debug
        try:
            adata = current_adata().copy()
            print(f"-> Working with adata: {adata.n_obs} cells, {adata.n_vars} genes")  # Debug
        
            # Basic QC calculations
            adata.var['mt'] = adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        
            print("-> QC metrics calculated")  # Debug
        
            qc_results.set({
                'adata': adata,
                'metrics': {
                    'n_cells_before': adata.n_obs,
                    'total_genes': adata.n_vars,
                    'mean_genes_per_cell': adata.obs['n_genes_by_counts'].mean()
                }
            })

            print(f"-> QC results set successfully")
            print(f"   - Cells: {adata.n_obs}")
            print(f"   - Genes: {adata.n_vars}")
            print(f"   - Mean genes/cell: {adata.obs['n_genes_by_counts'].mean():.2f}")
            print(f"QC results set: {qc_results.get()}")  # Debug
        
        except Exception as e:
            print(f"-> QC error: {e}")  # Debug
            import traceback
            traceback.print_exc()
        finally:
            is_processing.set(False)
            print("-> QC processing finished\n")  # Debug
    
    #@output
    @render.ui
    def qc_results_ui():
        print(">>> qc_results() called")
        results = qc_results.get()  # Get value from reactive.Value
        print(f">>> Has results: {results is not None}")
        
        if results is None:
            print(f">>> Returning 'Run QC' message.")
            return ui.p("Run QC analysis to see results.")
    
        metrics = results['metrics']
        print(f">>> Displaying metrics: Cells={metrics['n_cells_before']}")
        return ui.div(
            ui.h4("QC Metrics"),
            ui.p(f"-> Cells before QC: {metrics['n_cells_before']}"),
            ui.p(f"-> Total genes: {metrics['total_genes']}"),
            ui.p(f"-> Mean genes per cell: {metrics['mean_genes_per_cell']:.2f}"),
            #ui.output_plot("qc_plots")
        )

    @render.plot
    def qc_plots_module():
        print(">>> qc_plots_module() called")
        results = qc_results.get()
        print(f">>> Has results for plot: {results is not None}")

        if results is None:
            print(">>> Returning None - no results yet")
            return None
    
        adata = results['adata']
        print(f">>> Creating plot with {adata.n_obs} cells")
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(1, 3, figsize=(14, 5))

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
        if 'pct_counts_mt' in adata.obs.columns:
            axes[2].hist(adata.obs['pct_counts_mt'], bins=50, edgecolor='black')
            axes[2].set_xlabel('Mitochondrial %')
            axes[2].set_ylabel('Frequency')
            axes[2].set_title('Mitochondrial %')

        plt.tight_layout()
        print(f">>> Returning matplotlib figure ->")
        return fig

    return{
        'qc_results': qc_results
    }