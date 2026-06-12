from shiny import module, reactive, render, ui
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import muon as mu
import scanpy as sc
from sklearn.decomposition import NMF
import warnings
warnings.filterwarnings('ignore')

@module.ui
def create_multiomics_ui(module_id: str):
    return ui.div(
        ui.h3("Multi-omics Integration"),
        ui.layout_sidebar(
            ui.sidebar(
                ui.h5("Integration Parameters"),
                ui.input_select("integration_method", "Method",
                              choices=["WNN", "MOFA", "NMF", "DIABLO"],
                              selected="WNN"),
                ui.input_numeric("wnn_k", "WNN Neighbors", value=20, min=5, max=100),
                ui.input_numeric("factors", "Number of Factors", value=15, min=5, max=50),
                ui.input_select("modalities", "Modalities to Integrate", 
                              choices=[], multiple=True),
                ui.input_action_button("run_integration", "Run Integration", class_="btn-primary"),
                ui.hr(),
                ui.h6("Visualization"),
                ui.input_select("color_by", "Color By", choices=[]),
                ui.input_select("modality_plot", "Modality for Plot", choices=[]),
                width=350
            ),
            ui.div(
                ui.output_ui("integration_header"),
                ui.navset_pill(
                    ui.nav("Integration UMAP", ui.output_plot("integration_umap")),
                    ui.nav("Modality Weights", ui.output_plot("modality_weights")),
                    ui.nav("Factor Analysis", ui.output_plot("factor_analysis")),
                    ui.nav("Cross-modality", 
                          ui.input_select("gene_select", "Select Gene", choices=[]),
                          ui.input_select("peak_select", "Select Peak/Feature", choices=[]),
                          ui.output_plot("cross_modality_plot")
                    )
                ),
                class_="main-content"
            )
        )
    )

@module.server
def create_multiomics_server(module_id: str, current_adata, is_processing):
    integration_results = reactive.Value(None)
    mdata = reactive.Value(None)
    
    @reactive.Effect
    def _update_multiomics_ui():
        adata = current_adata()
        if adata is None:
            return
        
        # Check if we have MuData
        if hasattr(adata, 'mod') and isinstance(adata, mu.MuData):
            mdata.set(adata)
            modalities = list(adata.mod.keys())
            ui.update_select("modalities", choices=modalities)
            ui.update_select("modality_plot", choices=modalities)
            
            # Update color options
            all_obs_columns = set()
            for mod in modalities:
                if hasattr(adata.mod[mod], 'obs'):
                    all_obs_columns.update(adata.mod[mod].obs.columns)
            ui.update_select("color_by", choices=list(all_obs_columns))
    
    @reactive.Effect
    @reactive.event(input.run_integration)
    def _run_integration():
        if mdata() is None:
            return
        
        is_processing.set(True)
        try:
            current_mdata = mdata().copy()
            method = input.integration_method()
            modalities = input.modalities()
            
            if not modalities:
                raise ValueError("Please select at least one modality")
            
            if method == "WNN":
                # Weighted Nearest Neighbors integration
                mu.pp.intersect_obs(current_mdata)
                mu.tl.wnn(
                    current_mdata, 
                    k=input.wnk_k(),
                    modalities=modalities
                )
                
                # Compute UMAP on integrated space
                sc.tl.umap(current_mdata.obsm['X_wnn'])
                
                integration_results.set({
                    'method': 'WNN',
                    'mdata': current_mdata,
                    'modality_weights': current_mdata.uns['wnn']['weights'],
                    'embeddings': current_mdata.obsm['X_wnn']
                })
                
            elif method == "NMF":
                # Non-negative Matrix Factorization integration
                results = self._run_nmf_integration(current_mdata, modalities, input.factors())
                integration_results.set(results)
            
            # Update gene/feature selections
            self._update_feature_selections(current_mdata, modalities)
            
        except Exception as e:
            print(f"Integration error: {e}")
        finally:
            is_processing.set(False)
    
    def _run_nmf_integration(self, mdata, modalities, n_factors):
        """Run NMF-based integration"""
        # Concatenate features from selected modalities
        combined_features = []
        feature_names = []
        
        for modality in modalities:
            mod_data = mdata.mod[modality]
            if hasattr(mod_data, 'X'):
                # Use highly variable features if available
                if 'highly_variable' in mod_data.var:
                    features = mod_data[:, mod_data.var.highly_variable].X
                else:
                    # Use top expressed features
                    features = mod_data.X
                
                combined_features.append(features)
                feature_names.extend([f"{modality}_{gene}" for gene in mod_data.var_names])
        
        # Combine features
        from scipy import sparse
        if any(sparse.issparse(f) for f in combined_features):
            X_combined = sparse.hstack(combined_features)
        else:
            X_combined = np.hstack(combined_features)
        
        # Run NMF
        nmf = NMF(n_components=n_factors, random_state=42)
        W = nmf.fit_transform(X_combined)
        H = nmf.components_
        
        # Store results
        mdata.obsm['X_nmf'] = W
        mdata.varm['NMF_components'] = H.T
        
        return {
            'method': 'NMF',
            'mdata': mdata,
            'factors': n_factors,
            'reconstruction_error': nmf.reconstruction_err_,
            'feature_components': H
        }
    
    def _update_feature_selections(self, mdata, modalities):
        """Update gene and feature selection options"""
        gene_choices = []
        peak_choices = []
        
        for modality in modalities:
            mod_data = mdata.mod[modality]
            if 'rna' in modality.lower() or 'gex' in modality.lower():
                # RNA modality - treat as genes
                genes = mod_data.var_names.tolist()[:100]  # Top 100 genes
                gene_choices.extend([f"{modality}:{gene}" for gene in genes])
            elif 'atac' in modality.lower() or 'peak' in modality.lower():
                # ATAC modality - treat as peaks
                peaks = mod_data.var_names.tolist()[:100]  # Top 100 peaks
                peak_choices.extend([f"{modality}:{peak}" for peak in peaks])
        
        ui.update_select("gene_select", choices=gene_choices)
        ui.update_select("peak_select", choices=peak_choices)
    
    @render.ui
    def integration_header():
        if integration_results() is None:
            return ui.p("Select integration method and modalities to begin.")
        
        results = integration_results()
        method = results['method']
        
        stats_elements = [
            ui.value_box(
                title="Integration Method",
                value=method,
                theme="primary"
            )
        ]
        
        if 'modality_weights' in results:
            weights = results['modality_weights']
            for mod, weight in weights.items():
                stats_elements.append(
                    ui.value_box(
                        title=f"{mod} Weight",
                        value=f"{weight:.3f}",
                        theme="info"
                    )
                )
        
        return ui.div(
            ui.h4(f"Multi-omics Integration: {method}"),
            ui.layout_columns(*stats_elements)
        )
    
    @render.plot
    def integration_umap():
        if integration_results() is None:
            return None
        
        results = integration_results()
        mdata = results['mdata']
        color_by = input.color_by()
        
        # Get UMAP coordinates
        if 'X_umap' in mdata.obsm:
            umap_coords = mdata.obsm['X_umap']
        elif 'X_wnn' in mdata.obsm and 'X_umap' in mdata.obsm:
            umap_coords = mdata.obsm['X_umap']
        else:
            # Compute UMAP if not available
            sc.tl.umap(mdata)
            umap_coords = mdata.obsm['X_umap']
        
        # Get color values
        if color_by and color_by in mdata.obs:
            color_values = mdata.obs[color_by]
        else:
            color_values = None
        
        fig = px.scatter(
            x=umap_coords[:, 0],
            y=umap_coords[:, 1],
            color=color_values,
            title=f"Integrated UMAP - {results['method']}",
            labels={'x': 'UMAP1', 'y': 'UMAP2', 'color': color_by}
        )
        fig.update_layout(height=500)
        return fig
    
    @render.plot
    def modality_weights():
        if integration_results() is None or 'modality_weights' not in integration_results():
            return None
        
        weights = integration_results()['modality_weights']
        df = pd.DataFrame({
            'Modality': list(weights.keys()),
            'Weight': list(weights.values())
        })
        
        fig = px.bar(
            df, 
            x='Modality', 
            y='Weight',
            title="Modality Weights in Integration",
            color='Weight',
            color_continuous_scale='Viridis'
        )
        fig.update_layout(height=400)
        return fig
    
    @render.plot
    def cross_modality_plot():
        if integration_results() is None or not input.gene_select() or not input.peak_select():
            return None
        
        results = integration_results()
        mdata = results['mdata']
        
        gene = input.gene_select().split(":")[1]
        peak = input.peak_select().split(":")[1]
        gene_modality = input.gene_select().split(":")[0]
        peak_modality = input.peak_select().split(":")[0]
        
        # Get expression values
        gene_exp = mdata.mod[gene_modality][:, gene].X.flatten()
        peak_acc = mdata.mod[peak_modality][:, peak].X.flatten()
        
        # Create scatter plot
        fig = px.scatter(
            x=gene_exp,
            y=peak_acc,
            title=f"Cross-modality: {gene} vs {peak}",
            labels={'x': f'{gene} Expression', 'y': f'{peak} Accessibility'}
        )
        
        # Add correlation line
        correlation = np.corrcoef(gene_exp, peak_acc)[0, 1]
        fig.add_annotation(
            x=0.05, y=0.95,
            xref="paper", yref="paper",
            text=f"r = {correlation:.3f}",
            showarrow=False,
            bgcolor="white"
        )
        
        fig.update_layout(height=400)
        return fig