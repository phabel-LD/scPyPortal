import hashlib
import pickle
import asyncio
from pathlib import Path
from sqlmodel import SQLModel, Field, Session, select, create_engine
from typing import Optional, List, Dict, Any
from datetime import datetime, timedelta
import aiofiles
import psutil

from config import Config
from models import DatasetCache, PrecomputedResults, AnalysisSession

class CacheManager:
    """Advanced caching system with async support"""
    
    def __init__(self, database_url: str = Config.DATABASE_URL):
        self.engine = create_engine(database_url)
        SQLModel.metadata.create_all(self.engine)
    
    async def compute_file_hash(self, file_path: Path) -> str:
        """Compute SHA256 hash of file asynchronously"""
        sha256_hash = hashlib.sha256()
        async with aiofiles.open(file_path, "rb") as f:
            # Read file in chunks to handle large files
            while True:
                chunk = await f.read(8192)
                if not chunk:
                    break
                sha256_hash.update(chunk)
        return sha256_hash.hexdigest()
    
    async def cache_dataset(self, file_path: Path, original_filename: str) -> DatasetCache:
        """Cache a dataset with comprehensive metadata"""
        file_hash = await self.compute_file_hash(file_path)
        
        with Session(self.engine) as session:
            # Check if already cached
            existing = session.exec(
                select(DatasetCache).where(DatasetCache.file_hash == file_hash)
            ).first()
            
            if existing:
                return existing
            
            # Load AnnData to compute metrics
            import scanpy as sc
            adata = sc.read_h5ad(file_path)
            
            # Compute dataset metrics
            n_hvgs = None
            if 'highly_variable' in adata.var:
                n_hvgs = adata.var.highly_variable.sum()
            
            # Metadata summary
            metadata_summary = {
                'obs_columns': list(adata.obs.columns),
                'var_columns': list(adata.var.columns),
                'uns_keys': list(adata.uns.keys()),
                'obsm_keys': list(adata.obsm.keys()),
                'has_raw': adata.raw is not None
            }
            
            # Create cache entry
            cached_file_path = Config.DATASET_CACHE_DIR / f"{file_hash}.h5ad"
            adata.write_h5ad(cached_file_path)
            
            dataset_cache = DatasetCache(
                file_hash=file_hash,
                original_filename=original_filename,
                file_size_mb=file_path.stat().st_size / (1024 * 1024),
                n_cells=adata.n_obs,
                n_genes=adata.n_vars,
                n_hvgs=n_hvgs,
                pca_variance=[],
                metadata_summary=metadata_summary,
                processed_data_path=str(cached_file_path.relative_to(Config.BASE_DIR))
            )
            
            session.add(dataset_cache)
            session.commit()
            session.refresh(dataset_cache)
            
            return dataset_cache
    
    def get_cached_dataset(self, file_hash: str) -> Optional[DatasetCache]:
        """Retrieve cached dataset"""
        with Session(self.engine) as session:
            return session.exec(
                select(DatasetCache).where(DatasetCache.file_hash == file_hash)
            ).first()
    
    def store_precomputed_result(self, dataset_hash: str, computation_type: str, 
                               parameters: dict, results: dict):
        """Cache expensive computation results"""
        with Session(self.engine) as session:
            # Clean old results for same computation type
            old_results = session.exec(
                select(PrecomputedResults).where(
                    PrecomputedResults.dataset_hash == dataset_hash,
                    PrecomputedResults.computation_type == computation_type
                )
            ).all()
            
            for old in old_results:
                session.delete(old)
            
            # Store new results
            new_result = PrecomputedResults(
                dataset_hash=dataset_hash,
                computation_type=computation_type,
                parameters=parameters,
                results=results
            )
            
            session.add(new_result)
            session.commit()
    
    def get_precomputed_result(self, dataset_hash: str, computation_type: str, 
                             parameters: dict) -> Optional[dict]:
        """Retrieve cached computation results"""
        with Session(self.engine) as session:
            result = session.exec(
                select(PrecomputedResults).where(
                    PrecomputedResults.dataset_hash == dataset_hash,
                    PrecomputedResults.computation_type == computation_type
                )
            ).first()
            
            if result and result.parameters == parameters:
                return result.results
            
            return None
    
    def cleanup_old_cache(self, days_old: int = 30):
        """Remove cache entries older than specified days"""
        cutoff_date = datetime.now() - timedelta(days=days_old)
        
        with Session(self.engine) as session:
            # Find old datasets
            old_datasets = session.exec(
                select(DatasetCache).where(DatasetCache.upload_date < cutoff_date)
            ).all()
            
            for dataset in old_datasets:
                # Delete cached file
                cached_file = Config.BASE_DIR / dataset.processed_data_path
                if cached_file.exists():
                    cached_file.unlink()
                
                # Delete database entry
                session.delete(dataset)
            
            session.commit()