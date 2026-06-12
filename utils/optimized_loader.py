import numpy as np
import pandas as pd
import scanpy as sc
import h5py
import zarr
from pathlib import Path
import psutil
import warnings
from typing import Optional, Dict, Any
import threading
import asyncio

class OptimizedDataLoader:
    """Memory-efficient data loading and processing for large single-cell datasets"""
    
    def __init__(self, max_memory_gb: float = 8.0):
        self.max_memory_gb = max_memory_gb
        self.available_memory = self._get_available_memory()
        
    def _get_available_memory(self) -> float:
        """Get available memory in GB"""
        return psutil.virtual_memory().available / (1024 ** 3)
    
    def _check_memory_safety(self, expected_size_gb: float) -> bool:
        """Check if operation can be performed safely"""
        return (self.available_memory - expected_size_gb) > (self.max_memory_gb * 0.2)
    
    async def load_anndata_lazy(self, file_path: Path, **kwargs) -> sc.AnnData:
        """
        Load AnnData file with memory-efficient strategies
        """
        file_size_gb = file_path.stat().st_size / (1024 ** 3)
        
        # Strategy selection based on file size and available memory
        if file_size_gb > 2 or not self._check_memory_safety(file_size_gb * 3):
            return await self._load_large_dataset(file_path, **kwargs)
        else:
            return await self._load_standard(file_path, **kwargs)
    
    async def _load_standard(self, file_path: Path, **kwargs) -> sc.AnnData:
        """Standard loading for small to medium datasets"""
        try:
            adata = sc.read_h5ad(file_path, **kwargs)
            return adata
        except Exception as e:
            raise Exception(f"Failed to load {file_path}: {str(e)}")
    
    async def _load_large_dataset(self, file_path: Path, **kwargs) -> sc.AnnData:
        """Memory-efficient loading for large datasets"""
        
        # Strategy 1: Backed mode for very large datasets
        if file_path.stat().st_size > 500 * 1024 ** 2:  # > 500MB
            return await self._load_backed_mode(file_path, **kwargs)
        
        # Strategy 2: Chunked loading
        else:
            return await self._load_chunked(file_path, **kwargs)
    
    async def _load_backed_mode(self, file_path: Path, **kwargs) -> sc.AnnData:
        """Load in backed mode to keep data on disk"""
        try:
            # Load in backed mode
            adata = sc.read_h5ad(file_path, backed='r', **kwargs)
            
            # Create a lightweight in-memory version for operations
            # Only load metadata and small arrays
            lightweight_adata = sc.AnnData(
                obs=adata.obs.copy(),
                var=adata.var.copy(),
                uns=adata.uns.copy(),
                obsm=adata.obsm.copy() if hasattr(adata, 'obsm') else {},
                varm=adata.varm.copy() if hasattr(adata, 'varm') else {}
            )
            
            # Store reference to backed data
            lightweight_adata.uns['backed_file'] = str(file_path)
            lightweight_adata.uns['backed_mode'] = True
            
            return lightweight_adata
            
        except Exception as e:
            raise Exception(f"Backed mode loading failed: {str(e)}")
    
    async def _load_chunked(self, file_path: Path, chunk_size: int = 10000, **kwargs) -> sc.AnnData:
        """Load data in chunks for memory control"""
        try:
            with h5py.File(file_path, 'r') as f:
                # Read metadata first
                obs = pd.DataFrame({k: v[:] for k, v in f['obs'].items()})
                var = pd.DataFrame({k: v[:] for k, v in f['var'].items()})
                
                n_cells = obs.shape[0]
                n_genes = var.shape[0]
                
                # Initialize sparse matrix
                from scipy import sparse
                X = sparse.lil_matrix((n_cells, n_genes), dtype=np.float32)
                
                # Load in chunks
                for start_idx in range(0, n_cells, chunk_size):
                    end_idx = min(start_idx + chunk_size, n_cells)
                    
                    # Read chunk
                    chunk = f['X'][start_idx:end_idx]
                    if hasattr(chunk, 'shape') and len(chunk.shape) == 2:
                        X[start_idx:end_idx] = chunk
                    
                    # Yield control to event loop
                    await asyncio.sleep(0.001)
                
                # Convert to efficient format
                X = X.tocsr()
                
                adata = sc.AnnData(X=X, obs=obs, var=var)
                
                # Load other attributes
                if 'uns' in f:
                    adata.uns.update(self._load_dict(f['uns']))
                if 'obsm' in f:
                    for key in f['obsm'].keys():
                        adata.obsm[key] = f['obsm'][key][:]
                
                return adata
                
        except Exception as e:
            raise Exception(f"Chunked loading failed: {str(e)}")
    
    def _load_dict(self, h5_group):
        """Recursively load dictionaries from H5 groups"""
        result = {}
        for key in h5_group.keys():
            if isinstance(h5_group[key], h5py.Group):
                result[key] = self._load_dict(h5_group[key])
            else:
                result[key] = h5_group[key][:]
        return result
    
    def downsample_dataset(self, adata: sc.AnnData, target_cells: int = 50000, 
                          strategy: str = "random") -> sc.AnnData:
        """
        Downsample large datasets for responsive visualization
        """
        if adata.n_obs <= target_cells:
            return adata
        
        if strategy == "random":
            # Simple random sampling
            import random
            random.seed(42)
            indices = random.sample(range(adata.n_obs), target_cells)
            return adata[indices].copy()
        
        elif strategy == "cluster_aware":
            # Stratified sampling by clusters if available
            if 'leiden' in adata.obs:
                return self._cluster_aware_downsample(adata, target_cells)
            else:
                # Fallback to random
                return self.downsample_dataset(adata, target_cells, "random")
        
        else:
            raise ValueError(f"Unknown downsampling strategy: {strategy}")
    
    def _cluster_aware_downsample(self, adata: sc.AnnData, target_cells: int) -> sc.AnnData:
        """Downsample while preserving cluster proportions"""
        cluster_counts = adata.obs['leiden'].value_counts()
        sampled_indices = []
        
        for cluster, count in cluster_counts.items():
            cluster_indices = adata.obs[adata.obs['leiden'] == cluster].index
            n_sample = max(1, int(target_cells * (count / adata.n_obs)))
            
            if n_sample < len(cluster_indices):
                sampled = np.random.choice(cluster_indices, n_sample, replace=False)
                sampled_indices.extend(sampled)
            else:
                sampled_indices.extend(cluster_indices)
        
        # If we have too many cells, randomly sample to target
        if len(sampled_indices) > target_cells:
            sampled_indices = np.random.choice(sampled_indices, target_cells, replace=False)
        
        return adata[sampled_indices].copy()