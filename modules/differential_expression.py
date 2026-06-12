from shiny import module, reactive, render, ui
import pandas as pd
import numpy as np
import scanpy as sc
from statsmodels.stats.multitest import fdrcorrection
import warnings
warnings.filterwarnings('ignore')

@module.ui  # UNCOMMENT THIS
def create_de_ui():  # Remove module_id parameter
    return ui.div(
        ui.h3("Differential Expression Analysis"),
        ui.layout_sidebar(
            ui.sidebar(
                ui.h5("Analysis Parameters"),
                ui.input_select("group_by", "Group By", choices=[]),
                ui.input_select("group_1", "Group 1", choices=[], selected=None),
                ui.input_select("group_2", "Group 2", choices=[], selected=None),
                ui.input_select("test_method", "Statistical Test", 
                              choices=["t-test", "wilcoxon", "logreg"],
                              selected="wilcoxon"),
                ui.input_numeric("min_fold_change", "Min Log2 Fold Change", 
                               value=0.25, min=0, max=5, step=0.1),
                ui.input_numeric("max_pval", "Max Adjusted P-value", 
                               value=0.05, min=0, max=1, step=0.01),
                ui.input_numeric("min_cells", "Min Cells per Group", 
                               value=10, min=3, max=100),
                ui.input_action_button("run_de", "Run DE Analysis", class_="btn-primary"),
                ui.input_action_button("find_markers", "Find All Markers", class_="btn-success"),
                width=350
            ),
            ui.div(
                ui.output_ui("de_results_header"),
                ui.navset_pill(
                    ui.nav_panel("Volcano Plot", ui.output_plot("volcano_plot", height="500px")),
                    ui.nav_panel("DE Results", ui.output_data_frame("de_table")),
                    ui.nav_panel("Gene Expression", 
                          ui.input_select("de_gene_select", "Select Gene", choices=[]),
                          ui.output_plot("gene_expression_plot", height="400px")
                    ),
                    ui.nav_panel("All Markers",
                                 ui.input_select("marker_group_select", "Select Cell Type", choices=[]),
                                 ui.output_data_frame("all_markers_table")
                    )
                ),
                class_="main-content"
            )
        )
    )

@module.server
def create_de_server(input, output, session, current_adata, is_processing):
    de_results = reactive.Value(None)
    marker_genes = reactive.Value(None)
    
    # Update UI when data changes
    @reactive.Effect
    def _update_de_ui():
        adata = current_adata()
        if adata is None:
            return
        
        # Update group_by choices
        categorical_cols = [col for col in adata.obs.columns 
                          if adata.obs[col].dtype.name == 'category' or 
                          adata.obs[col].nunique() < 20]
        ui.update_select("group_by", choices=categorical_cols)
    
    @reactive.Effect
    @reactive.event(input.group_by)
    def _update_group_choices():
        adata = current_adata()
        group_by = input.group_by()
        if adata is None or not group_by:
            return
        
        groups = sorted(adata.obs[group_by].unique().astype(str))
        ui.update_select("group_1", choices=groups)
        ui.update_select("group_2", choices=groups)
    
    @reactive.Effect
    @reactive.event(input.run_de)
    def _run_differential_expression():
        print("\n" + "="*60)
        print("DE ANALYSIS BUTTON CLICKED")
        print("="*60)
        
        if current_adata() is None:
            print("❌ No data loaded")
            return
        
        is_processing.set(True)
        print("🔄 Starting DE analysis...")
        try:
            adata = current_adata().copy()
            group_by = input.group_by()
            group_1 = input.group_1()
            group_2 = input.group_2()
            
            print(f"  Comparing: {group_1} vs {group_2} (by {group_by})")
            
            if not all([group_by, group_1, group_2]):
                raise ValueError("Please select both groups for comparison")
            
            # Filter to groups of interest
            mask = adata.obs[group_by].astype(str).isin([group_1, group_2])
            adata_subset = adata[mask].copy()
            
            print(f"  Subset: {adata_subset.n_obs} cells")
            
            # Ensure we have enough cells
            group_counts = adata_subset.obs[group_by].value_counts()
            if any(group_counts < input.min_cells()):
                raise ValueError(f"Groups must have at least {input.min_cells()} cells")
            
            # Run differential expression
            print("  → Running differential expression test...")
            sc.tl.rank_genes_groups(
                adata_subset, 
                groupby=group_by, 
                groups=[group_1],
                reference=group_2,
                method=input.test_method(),
                use_raw=False
            )
            
            # Extract results
            de_df = sc.get.rank_genes_groups_df(adata_subset, group=group_1)
            
            # Apply multiple testing correction
            de_df['pvals_adj'] = fdrcorrection(de_df['pvals'])[1]
            
            # Filter results
            de_df = de_df[
                (de_df['pvals_adj'] <= input.max_pval()) & 
                (abs(de_df['logfoldchanges']) >= input.min_fold_change())
            ].sort_values('scores', ascending=False)
            
            print(f"✅ Found {len(de_df)} significant genes")
            
            de_results.set({
                'data': de_df,
                'parameters': {
                    'group_by': group_by,
                    'group_1': group_1,
                    'group_2': group_2,
                    'test_method': input.test_method()
                },
                'adata_subset': adata_subset
            })
            
            # Update gene selection
            top_genes = de_df.head(50)['names'].tolist()
            ui.update_select("de_gene_select", choices=top_genes)
            
        except Exception as e:
            print(f"❌ DE analysis error: {e}")
            import traceback
            traceback.print_exc()
        finally:
            is_processing.set(False)
            print("✅ DE analysis finished\n")
    
    @reactive.Effect
    @reactive.event(input.find_markers)
    def _find_all_markers():
        print("\n" + "="*60)
        print("FIND ALL MARKERS BUTTON CLICKED")
        print("="*60)
        
        if current_adata() is None:
            print("❌ No data loaded")
            return
        
        is_processing.set(True)
        print("🔄 Finding markers for all groups...")
        try:
            adata = current_adata().copy()
            group_by = input.group_by()
            
            if not group_by:
                raise ValueError("Please select a grouping variable")
            
            # Find markers for all groups
            sc.tl.rank_genes_groups(
                adata, 
                groupby=group_by, 
                method=input.test_method(),
                use_raw=False
            )
            
            # Extract results for all groups
            marker_dict = {}
            for group in adata.obs[group_by].cat.categories:
                de_df = sc.get.rank_genes_groups_df(adata, group=group)
                de_df['pvals_adj'] = fdrcorrection(de_df['pvals'])[1]
                
                # Filter and store
                marker_dict[group] = de_df[
                    (de_df['pvals_adj'] <= input.max_pval()) & 
                    (abs(de_df['logfoldchanges']) >= input.min_fold_change())
                ].head(100)  # Top 100 markers per group
                
                print(f"  {group}: {len(marker_dict[group])} markers")
            
            marker_genes.set(marker_dict)
            print(f"✅ Marker finding complete!")
            
        except Exception as e:
            print(f"❌ Marker finding error: {e}")
            import traceback
            traceback.print_exc()
        finally:
            is_processing.set(False)
            print("✅ Marker finding finished\n")
    
    @render.ui
    def de_results_header():
        results = de_results.get()
        if results is None:
            return ui.p("Configure parameters and run analysis to see results.")
        
        params = results['parameters']
        n_genes = len(results['data'])
        
        return ui.div(
            ui.h4(f"Differential Expression: {params['group_1']} vs {params['group_2']}"),
            ui.layout_columns(
                ui.value_box(
                    title="Significant Genes",
                    value=n_genes,
                    theme="primary"
                ),
                ui.value_box(
                    title="Test Method",
                    value=params['test_method'],
                    theme="info"
                ),
                ui.value_box(
                    title="Comparison",
                    value=f"{params['group_1']} vs {params['group_2']}",
                    theme="success"
                )
            )
        )
    
    @render.plot
    def volcano_plot():
        print(">>> volcano_plot() called")
        results = de_results.get()
        if results is None:
            return None
        
        de_df = results['data'].copy()
        
        import matplotlib.pyplot as plt
        
        # Create volcano plot data
        de_df['-log10(padj)'] = -np.log10(de_df['pvals_adj'] + 1e-300)  # Avoid log(0)
        de_df['significant'] = (de_df['pvals_adj'] <= 0.05) & (abs(de_df['logfoldchanges']) >= 0.5)
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Plot non-significant in gray
        non_sig = de_df[~de_df['significant']]
        ax.scatter(non_sig['logfoldchanges'], non_sig['-log10(padj)'], 
                  c='gray', s=10, alpha=0.5, label='Not significant')
        
        # Plot significant in red
        sig = de_df[de_df['significant']]
        ax.scatter(sig['logfoldchanges'], sig['-log10(padj)'], 
                  c='red', s=15, alpha=0.7, label='Significant')
        
        # Add threshold lines
        ax.axhline(y=-np.log10(0.05), linestyle='--', color='black', linewidth=1, alpha=0.5)
        ax.axvline(x=0.5, linestyle='--', color='black', linewidth=1, alpha=0.5)
        ax.axvline(x=-0.5, linestyle='--', color='black', linewidth=1, alpha=0.5)
        
        ax.set_xlabel('Log2 Fold Change', fontsize=10)
        ax.set_ylabel('-Log10(Adj. P-value)', fontsize=10)
        ax.set_title('Volcano Plot', fontsize=12)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        print(">>> Returning volcano plot ✅")
        return fig
    
    @render.data_frame
    def de_table():
        results = de_results.get()
        if results is None:
            return render.DataGrid(pd.DataFrame(), height=400)
        
        de_df = results['data'].copy()
        # Format for display
        display_df = de_df.head(1000)[['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'scores']]
        display_df.columns = ['Gene', 'Log2FC', 'P-value', 'Adj. P-value', 'Score']
        display_df = display_df.round(4)
        
        return render.DataGrid(display_df, height=400, filters=True)
    
    @render.plot
    def gene_expression_plot():
        print(">>> gene_expression_plot() called")
        results = de_results.get()
        gene = input.de_gene_select()
        
        if results is None or not gene:
            return None
        
        adata = results['adata_subset']
        params = results['parameters']
        
        import matplotlib.pyplot as plt
        
        # Get expression data
        gene_expr = adata[:, gene].X
        if hasattr(gene_expr, 'toarray'):
            gene_expr = gene_expr.toarray().flatten()
        else:
            gene_expr = gene_expr.flatten()
        
        # Create violin plot
        fig, ax = plt.subplots(figsize=(6, 5))
        
        groups = adata.obs[params['group_by']].unique()
        positions = []
        data_to_plot = []
        
        for i, group in enumerate(sorted(groups)):
            mask = adata.obs[params['group_by']] == group
            data_to_plot.append(gene_expr[mask])
            positions.append(i)
        
        parts = ax.violinplot(data_to_plot, positions=positions, showmeans=True, showmedians=True)
        
        ax.set_xticks(positions)
        ax.set_xticklabels(sorted(groups))
        ax.set_ylabel('Expression', fontsize=10)
        ax.set_xlabel(params['group_by'], fontsize=10)
        ax.set_title(f'Expression of {gene}', fontsize=12)
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        print(">>> Returning gene expression plot ✅")
        return fig
    
    @reactive.Effect
    def _update_marker_group_options():
        """Update group options after finding markers"""
        markers = marker_genes.get()
        if markers is None:
            return
    
        ui.update_select("marker_group_select", choices=list(markers.keys()))

    @render.data_frame
    def all_markers_table():
        markers = marker_genes.get()
        group = input.marker_group_select()
    
        if markers is None or not group or group not in markers:
            return render.DataGrid(pd.DataFrame(), height=400)
    
        marker_df = markers[group].copy()
        display_df = marker_df[['names', 'logfoldchanges', 'pvals_adj', 'scores']].head(100)
        display_df.columns = ['Gene', 'Log2FC', 'Adj. P-value', 'Score']
        display_df = display_df.round(4)
    
        return render.DataGrid(display_df, height=500, filters=True)

    return {
        'de_results': de_results,
        'marker_genes': marker_genes
    }