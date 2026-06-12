from shiny import module, reactive, render, ui
import pandas as pd
import numpy as np
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

@module.ui  # UNCOMMENT THIS
def create_trajectory_ui():  # Remove module_id parameter
    return ui.div(
        ui.h3("Trajectory Inference & Pseudotime Analysis"),
        ui.layout_sidebar(
            ui.sidebar(
                ui.h5("Trajectory Parameters"),
                ui.input_select("trajectory_method", "Method", 
                              choices=["diffusion_map", "dpt"],
                              selected="dpt"),
                ui.input_select("root_grouping", "Root Cell Grouping", choices=[]),
                ui.input_select("root_group_value", "Root Cluster/Type", choices=[]),
                ui.input_select("color_groups", "Color By", choices=[]),
                ui.input_numeric("trajectory_n_dcs", "Diffusion Components", value=10, min=3, max=50),
                ui.input_numeric("trajectory_n_neighbors", "Number of Neighbors", value=30, min=5, max=100),
                ui.input_action_button("run_trajectory", "Run Trajectory Analysis", class_="btn-primary"),
                width=350
            ),
            ui.div(
                ui.output_ui("trajectory_header"),
                ui.navset_pill(
                    ui.nav_panel("Diffusion Map", ui.output_plot("diffusion_plot", height="500px")),
                    ui.nav_panel("Pseudotime", ui.output_plot("pseudotime_plot", height="500px")),
                    ui.nav_panel("Gene Trends", 
                          ui.input_select("trend_gene", "Select Gene", choices=[]),
                          ui.output_plot("gene_trend_plot", height="400px")
                    )
                ),
                class_="main-content"
            )
        )
    )

@module.server
def create_trajectory_server(input, output, session, current_adata, is_processing):
    trajectory_results = reactive.Value(None)
    trajectory_adata = reactive.Value(None)
    
    @reactive.Effect
    def _update_trajectory_ui():
        """Update UI options when data is available"""
        adata = current_adata()
        if adata is None:
            print(">>> No data available yet")
            return
        
        print(">>> Updating trajectory UI options")

        # Update group choices for root cell selection
        categorical_cols = [col for col in adata.obs.columns 
                          if adata.obs[col].dtype.name == 'category' or 
                          adata.obs[col].nunique() < 20]
        
        print(f">>> Found categorical columns: {categorical_cols}")
        
        if len(categorical_cols) > 0:
            ui.update_select("root_grouping", choices=categorical_cols)
            ui.update_select("color_groups", choices=categorical_cols)
            print(f">>> Updated dropdowns with {len(categorical_cols)} options")

        # Update gene choices
        if 'highly_variable' in adata.var:
            gene_choices = adata.var_names[adata.var.highly_variable].tolist()[:200]
        else:
            gene_choices = adata.var_names.tolist()[:200]
        ui.update_select("trend_gene", choices=gene_choices)
    
        print(f">>> UI updated with {len(categorical_cols)} grouping options")

    @reactive.Effect
    @reactive.event(input.root_grouping)
    def _update_root_group_values():
        """Update root group value choices when grouping changes"""
        print(f">>> root_grouping changed to: {input.root_grouping()}")
        
        adata = current_adata()
        grouping = input.root_grouping()
    
        if adata is None or not grouping:
            print(">>> No data or no grouping selected")
            return
    
        # Get unique values in the selected grouping
        values = sorted(adata.obs[grouping].unique().astype(str))
        print(f">>> Found {len(values)} unique values in {grouping}: {values}")
        ui.update_select("root_group_value", choices=values)

    @reactive.Effect
    @reactive.event(input.run_trajectory)
    def _run_trajectory_analysis():
        print("\n" + "="*60)
        print("TRAJECTORY ANALYSIS BUTTON CLICKED")
        print("="*60)
        
        if current_adata() is None:
            print("❌ No data loaded")
            return
        
        adata = current_adata()
        categorical_cols = [col for col in adata.obs.columns
                            if adata.obs[col].dtype.name == 'category' or
                            adata.obs[col].nunique() < 20]
        
        if len(categorical_cols) > 0:
            print(f"Forcing UI update with columns: {categorical_cols}")
            ui.update_select("root_grouping", choices=categorical_cols)
            ui.update_select("color_groups", choices=categorical_cols)

        is_processing.set(True)
        print("🔄 Starting trajectory analysis...")
        try:
            adata = current_adata().copy()
            method = input.trajectory_method()
            
            print(f"  Method: {method}")
            print(f"  Processing {adata.n_obs} cells")
            
            # Ensure neighbors are computed
            if 'neighbors' not in adata.uns:
                print("  → Computing neighbors...")
                sc.pp.neighbors(adata, n_neighbors=input.trajectory_n_neighbors())
            
            if method == "diffusion_map":
                print("  → Computing diffusion map...")
                sc.tl.diffmap(adata, n_comps=input.trajectory_n_dcs())
                
                trajectory_results.set({
                    'method': 'diffusion_map',
                    'components': adata.obsm['X_diffmap'],
                    'eigenvalues': adata.uns['diffmap_evals']
                })
                
                print("✅ Diffusion map complete!")
                
            elif method == "dpt":
                print("  → Computing DPT (diffusion pseudotime)...")
                
                # Get root cell from selected cluster
                root_grouping = input.root_grouping()
                root_value = input.root_group_value()
                
                if not root_grouping or not root_value:
                    raise ValueError("Please select both a grouping variable and a specific cluster/type for the root cell")
            
                # Find all cells in the selected group
                root_mask = adata.obs[root_grouping].astype(str) == root_value
                root_cell_indices = np.where(root_mask)[0]

                if len(root_cell_indices) == 0:
                    raise ValueError(f"No cells found in {root_grouping}={root_value}")
                
                # Randomly select one cell from this group as root
                root_cell_index = np.random.choice(root_cell_indices)
                
                print(f"  Root group: {root_grouping}")
                print(f"  Root cell index: {root_cell_index}")
                print(f"  Available cells in group: {len(root_cell_indices)}")
                print(f"  Selected root cell index: {root_cell_index}")
                
                adata.uns['iroot'] = root_cell_index
                
                sc.tl.diffmap(adata, n_comps=input.trajectory_n_dcs())
                sc.tl.dpt(adata)
                
                trajectory_results.set({
                    'method': 'dpt',
                    'pseudotime': adata.obs['dpt_pseudotime'],
                    'root_cell': root_cell_index,
                    'root_group': root_value,
                    'components': adata.obsm['X_diffmap']
                })
                
                print("✅ DPT complete!")
            
            trajectory_adata.set(adata)
            
        except Exception as e:
            print(f"❌ Trajectory analysis error: {e}")
            import traceback
            traceback.print_exc()
        finally:
            is_processing.set(False)
            print("✅ Trajectory analysis finished\n")
    
    @render.ui
    def trajectory_header():
        results = trajectory_results.get()
        if results is None:
            return ui.p("Configure parameters and run trajectory analysis.")
        
        method = results['method']
        
        return ui.div(
            ui.h4(f"Trajectory Analysis: {method.upper()}"),
            ui.layout_columns(
                ui.value_box(
                    title="Method",
                    value=method.upper(),
                    theme="primary"
                ),
                ui.value_box(
                    title="Status",
                    value="Completed",
                    theme="success"
                )
            )
        )
    
    @render.plot
    def diffusion_plot():
        print(">>> diffusion_plot() called")
        results = trajectory_results.get()
        adata = trajectory_adata.get()
        
        if results is None or adata is None:
            return None
        
        components = results['components']
        color_by = input.color_groups()
        
        print(f">>> Creating diffusion map, colored by: {color_by}")
        
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        if color_by and color_by in adata.obs.columns:
            # Color by selected grouping
            for category in sorted(adata.obs[color_by].unique()):
                mask = adata.obs[color_by] == category
                ax.scatter(
                    components[mask, 0], 
                    components[mask, 1],
                    label=str(category),
                    s=15,
                    alpha=0.7
                )
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        else:
            # No coloring
            ax.scatter(components[:, 0], components[:, 1], s=15, alpha=0.7)
        
        ax.set_xlabel('DC1', fontsize=10)
        ax.set_ylabel('DC2', fontsize=10)
        ax.set_title('Diffusion Map Components 1-2', fontsize=12)
        
        plt.tight_layout()
        print(">>> Returning diffusion map ✅")
        return fig
    
    @render.plot
    def pseudotime_plot():
        print(">>> pseudotime_plot() called")
        results = trajectory_results.get()
        adata = trajectory_adata.get()
        
        if results is None or adata is None or results['method'] != 'dpt':
            return None
        
        print(">>> Creating pseudotime plot")
        
        import matplotlib.pyplot as plt
        
        # Use UMAP for visualization if available, otherwise use diffusion map
        if 'X_umap' in adata.obsm:
            coords = adata.obsm['X_umap']
            xlabel, ylabel = 'UMAP1', 'UMAP2'
        else:
            coords = results['components']
            xlabel, ylabel = 'DC1', 'DC2'
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        scatter = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=results['pseudotime'],
            cmap='viridis',
            s=20,
            alpha=0.8
        )
        
        plt.colorbar(scatter, ax=ax, label='Pseudotime')
        
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title('Pseudotime Trajectory', fontsize=12)
        
        plt.tight_layout()
        print(">>> Returning pseudotime plot ✅")
        return fig
    
    @render.plot
    def gene_trend_plot():
        print(">>> gene_trend_plot() called")
        results = trajectory_results.get()
        adata = trajectory_adata.get()
        gene = input.trend_gene()
        
        if results is None or adata is None or results['method'] != 'dpt' or not gene:
            return None
        
        print(f">>> Creating gene trend plot for {gene}")
        
        import matplotlib.pyplot as plt
        
        # Get gene expression
        if gene not in adata.var_names:
            print(f"⚠️ Gene {gene} not found")
            return None
        
        gene_expr = adata[:, gene].X
        if hasattr(gene_expr, 'toarray'):
            gene_expr = gene_expr.toarray().flatten()
        else:
            gene_expr = gene_expr.flatten()
        
        # Plot gene expression along pseudotime
        pseudotime = results['pseudotime']
        
        # Create trend plot with binning
        fig, ax = plt.subplots(figsize=(8, 5))
        
        # Scatter plot
        ax.scatter(pseudotime, gene_expr, alpha=0.3, s=10, color='gray')
        
        # Add smoothed trend line
        from scipy.ndimage import uniform_filter1d
        sorted_idx = np.argsort(pseudotime)
        smoothed = uniform_filter1d(gene_expr[sorted_idx], size=50)
        ax.plot(pseudotime[sorted_idx], smoothed, color='red', linewidth=2, label='Trend')
        
        ax.set_xlabel('Pseudotime', fontsize=10)
        ax.set_ylabel('Expression', fontsize=10)
        ax.set_title(f'Expression of {gene} along Pseudotime', fontsize=12)
        ax.legend()
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        print(">>> Returning gene trend plot ✅")
        return fig
    
    return {
        'trajectory_results': trajectory_results
    }