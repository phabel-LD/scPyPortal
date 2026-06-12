from pathlib import Path
from typing import Dict, Any
import os

class Config:
    """Application configuration"""
    
    # Paths
    BASE_DIR = Path(__file__).parent
    CACHE_DIR = BASE_DIR / "cache"
    DATASET_CACHE_DIR = CACHE_DIR / "datasets"
    TEMP_DIR = CACHE_DIR / "temp"
    
    # Database
    DATABASE_URL = "sqlite:///./scpyportal.db"
    
    # Analysis defaults
    DEFAULT_N_NEIGHBORS = 15
    DEFAULT_RESOLUTION = 0.5
    DEFAULT_N_PCS = 50
    DEFAULT_PERPLEXITY = 30
    
    # Performance
    MAX_FILE_SIZE_MB = 500
    CACHE_EXPIRY_DAYS = 30
    BACKGROUND_PROCESSING = True
    
    # Visualization
    PLOTLY_THEME = "plotly_white"
    COLOR_SCHEMES = ["viridis", "plasma", "inferno", "magma", "cividis"]
    
    @classmethod
    def setup_directories(cls):
        """Create necessary directories"""
        cls.DATASET_CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cls.TEMP_DIR.mkdir(parents=True, exist_ok=True)

# Initialize directories
Config.setup_directories()