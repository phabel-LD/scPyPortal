import pandas as pd
import numpy as np
from pathlib import Path
import scanpy as sc

def validate_anndata_file(file_path: Path) -> dict:
    """
    Validate an AnnData file before processing
    Returns: dict with 'is_valid' and 'message'
    """
    try:
        # Quick check without full load
        with sc.read_h5ad(file_path, backed='r') as adata:
            validation_result = {
                'is_valid': True,
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'has_X': adata.X is not None,
                'obs_columns': list(adata.obs.columns) if hasattr(adata, 'obs') else [],
                'var_columns': list(adata.var.columns) if hasattr(adata, 'var') else []
            }
            
            # Additional checks
            if adata.n_obs == 0 or adata.n_vars == 0:
                validation_result['is_valid'] = False
                validation_result['message'] = "Dataset has zero cells or genes"
            elif adata.n_obs > 1000000:
                validation_result['message'] = "Warning: Very large dataset (>1M cells)"
            else:
                validation_result['message'] = "Dataset looks good"
                
        return validation_result
        
    except Exception as e:
        return {
            'is_valid': False,
            'message': f"Invalid AnnData file: {str(e)}"
        }

def validate_parameters(parameters: dict) -> list:
    """Validate analysis parameters"""
    errors = []
    
    if 'n_neighbors' in parameters and parameters['n_neighbors'] <= 0:
        errors.append("Number of neighbors must be positive")
    
    if 'resolution' in parameters and parameters['resolution'] <= 0:
        errors.append("Clustering resolution must be positive")
    
    if 'n_pcs' in parameters and (parameters['n_pcs'] < 1 or parameters['n_pcs'] > 100):
        errors.append("Number of PCs must be between 1 and 100")
    
    return errors