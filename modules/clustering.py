from shiny import module, reactive, render, ui
import scanpy as sc
import numpy as np
from sklearn.metrics import silhouette_score

@module.ui  # UNCOMMENT THIS
def create_clustering_ui():  # Remove module_id parameter
    return ui.div(
        ui.h3("Dimensionality Reduction & Clustering"),
        ui.layout_sidebar(
            ui.sidebar(
                ui.h5("Parameters"),
                ui.input_numeric("n_pcs", "Number of PCs", value=50, min=10, max=100),
                ui.input_numeric("clustering_n_neighbors", "Number of Neighbors", value=15, min=5, max=100),
                ui.input_numeric("resolution", "Clustering Resolution", value=0.5, min=0.1, max=2.0, step=0.1),
                ui.input_select("metric", "Distance Metric", 
                              choices=["euclidean", "manhattan", "cosine"],
                              selected="euclidean"),
                ui.input_numeric("random_state", "Random Seed", value=42),
                ui.input_action_button("run_clustering", "Run Clustering", class_="btn-primary"),
                ui.hr(),
                ui.h5("Visualization Options"),
                ui.input_select("color_by", "Color UMAP by:",
                                choices=["leiden", "CellType", "Cell Type"],
                                selected="leiden"),
                width=300
            ),
            ui.div(
                ui.output_ui("clustering_results_ui"),  # Renamed
                ui.output_plot("umap_plot", height="500px"),
                ui.output_plot("cluster_sizes_plot", height="400px"),
                class_="main-content"
            )
        )
    )

@module.server
def create_clustering_server(input, output, session, current_adata, is_processing):
    adata_processed = reactive.Value(None)
    clustering_results = reactive.Value(None)
    
    @reactive.Effect
    @reactive.event(input.run_clustering)
    def _run_clustering():
        print("\n" + "="*60)
        print("CLUSTERING BUTTON CLICKED")
        print("="*60)
        
        if current_adata() is None:
            print("❌ No data loaded")
            return
        
        is_processing.set(True)
        print("🔄 Starting clustering...")
        try:
            adata = current_adata().copy()
            print(f"✅ Processing {adata.n_obs} cells, {adata.n_vars} genes")
            
            # Run clustering pipeline
            print("  → Running PCA...")
            sc.tl.pca(adata, n_comps=input.n_pcs(), random_state=input.random_state())
            
            print("  → Computing neighbors...")
            sc.pp.neighbors(
                adata, 
                n_neighbors=input.clustering_n_neighbors(),
                metric=input.metric(),
                random_state=input.random_state()
            )
            
            print("  → Computing UMAP...")
            sc.tl.umap(adata, random_state=input.random_state())
            
            print("  → Running Leiden clustering...")
            sc.tl.leiden(adata, resolution=input.resolution(), random_state=input.random_state())
            
            # Calculate clustering metrics
            print("  → Calculating silhouette score...")
            silhouette_avg = silhouette_score(
                adata.obsm['X_pca'][:, :input.n_pcs()],
                adata.obs['leiden'].astype('category').cat.codes
            )
            
            adata_processed.set(adata)
            clustering_results.set({
                'silhouette_score': silhouette_avg,
                'n_clusters': len(adata.obs['leiden'].unique()),
                'parameters': {
                    'n_pcs': input.n_pcs(),
                    'n_neighbors': input.clustering_n_neighbors(),
                    'resolution': input.resolution()
                }
            })
            
            print(f"✅ Clustering complete! Found {len(adata.obs['leiden'].unique())} clusters")
            print(f"   Silhouette score: {silhouette_avg:.3f}")
            
        except Exception as e:
            print(f"❌ Clustering error: {e}")
            import traceback
            traceback.print_exc()
        finally:
            is_processing.set(False)
            print("✅ Clustering finished\n")
    
    @render.ui
    def clustering_results_ui():  # Renamed to avoid collision
        print(">>> clustering_results_ui() called")
        
        if adata_processed.get() is None:
            return ui.p("Click 'Run Clustering' to analyze your data.")
        
        adata = adata_processed.get()
        results = clustering_results.get()
        
        return ui.div(
            ui.h4("Clustering Results (Leiden Algorithm)"),
            ui.p("Note: Metrics below are for the computed Leiden clusters, not the visualization coloring.",
                 style="font-size: 0.9em; color: #666; font-style: italic;"),
            ui.layout_columns(
                ui.value_box(
                    title="Leiden Clusters Found",
                    value=results['n_clusters'],
                    theme="primary"
                ),
                ui.value_box(
                    title="Silhouette Score (Leiden)",
                    value=f"{results['silhouette_score']:.3f}",
                    theme="success"
                )
            ),
        )
    
    @render.plot
    def umap_plot():
        print(">>> umap_plot() called")
        if adata_processed.get() is None:
            return None
        
        adata = adata_processed.get()
        color_by = input.color_by()

        # Check if the selected column exists
        if color_by not in adata.obs.columns:
            print(f"⚠️ Column '{color_by}' not found, using 'leiden'")
            color_by = 'leiden'

        print(f">>> Creating UMAP plot with {len(adata.obs['leiden'].unique())} clusters")
        
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Get unique categories
        categories = sorted(adata.obs[color_by].unique())

        # Plot UMAP colored by selected variable
        for category  in categories:
            mask = adata.obs[color_by] == category
            ax.scatter(
                adata.obsm['X_umap'][mask, 0],
                adata.obsm['X_umap'][mask, 1],
                label=f'Cluster {category}',
                s=15,
                alpha=0.7
            )
        
        ax.set_xlabel('UMAP 1', fontsize=10)
        ax.set_ylabel('UMAP 2', fontsize=10)
        ax.set_title(f'UMAP Projection - Colored by {color_by}', fontsize=12)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        
        plt.tight_layout()
        print(">>> Returning UMAP figure ✅")
        return fig
    
    @render.plot
    def cluster_sizes_plot():
        print(">>> cluster_sizes_plot() called")
        if adata_processed.get() is None:
            return None
        
        adata = adata_processed.get()
        
        import matplotlib.pyplot as plt
        
        # Count cells per cluster
        cluster_counts = adata.obs['leiden'].value_counts().sort_index()
        
        fig, ax = plt.subplots(figsize=(10, 5))
        cluster_counts.plot(kind='bar', ax=ax, color='steelblue', edgecolor='black')
        ax.set_xlabel('Cluster')
        ax.set_ylabel('Number of Cells')
        ax.set_title('Cells per Cluster')
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        print(">>> Returning cluster sizes figure ✅")
        return fig
    
    return {
        'adata_processed': adata_processed,
        'clustering_results': clustering_results
    }