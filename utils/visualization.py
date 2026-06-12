import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots

def create_umap_plot(adata, color_by=None, title="UMAP Projection"):
    """Create UMAP visualization"""
    if 'X_umap' not in adata.obsm:
        raise ValueError("UMAP coordinates not found")
    
    umap_df = pd.DataFrame({
        'UMAP1': adata.obsm['X_umap'][:, 0],
        'UMAP2': adata.obsm['X_umap'][:, 1]
    })
    
    if color_by and color_by in adata.obs:
        umap_df['color'] = adata.obs[color_by]
        fig = px.scatter(umap_df, x='UMAP1', y='UMAP2', color='color', title=title)
    else:
        fig = px.scatter(umap_df, x='UMAP1', y='UMAP2', title=title)
    
    fig.update_layout(height=500, showlegend=True)
    return fig

def create_violin_plot(adata, gene, group_by=None, title="Gene Expression"):
    """Create violin plot for gene expression"""
    if gene not in adata.var_names:
        raise ValueError(f"Gene {gene} not found in dataset")
    
    expression = adata[:, gene].X.flatten() if hasattr(adata[:, gene].X, 'flatten') else adata[:, gene].X.toarray().flatten()
    
    plot_df = pd.DataFrame({'expression': expression})
    
    if group_by and group_by in adata.obs:
        plot_df['group'] = adata.obs[group_by]
        fig = px.violin(plot_df, x='group', y='expression', box=True, title=title)
    else:
        fig = px.violin(plot_df, y='expression', box=True, title=title)
    
    fig.update_layout(height=400)
    return fig

def create_heatmap_plot(adata, genes, group_by, title="Expression Heatmap"):
    """Create a heatmap for multiple genes"""
    # Subset to selected genes
    valid_genes = [g for g in genes if g in adata.var_names]
    if not valid_genes:
        raise ValueError("No valid genes found")
    
    # Get expression matrix
    expr_matrix = adata[:, valid_genes].X
    if hasattr(expr_matrix, 'toarray'):
        expr_matrix = expr_matrix.toarray()
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=expr_matrix,
        x=valid_genes,
        y=adata.obs[group_by] if group_by in adata.obs else list(range(adata.n_obs)),
        colorscale='Viridis'
    ))
    
    fig.update_layout(
        title=title,
        height=600,
        xaxis_title="Genes",
        yaxis_title="Cells" if group_by not in adata.obs else group_by
    )
    
    return fig