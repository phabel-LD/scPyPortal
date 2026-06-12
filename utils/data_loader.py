from pathlib import Path
import scanpy as sc
from .cache_manager import CacheManager
from config import Config

class DataLoader:
    def __init__(self, cache_manager: CacheManager):
        self.cache_manager = cache_manager
    
    def load_cached_dataset(self, file_hash: str):
        """Load dataset from cache"""
        cached_dataset = self.cache_manager.get_cached_dataset(file_hash)
        if cached_dataset and Path(cached_dataset.processed_data_path).exists():
            return sc.read_h5ad(cached_dataset.processed_data_path)
        return None
    
    async def load_dataset(self, file_info):
        """Main method to load dataset (used in app.py)"""
        file_path = Path(file_info["datapath"])
        dataset_cache = await self.cache_manager.cache_dataset(file_path, file_info["name"])
        return self.load_cached_dataset(dataset_cache.file_hash)