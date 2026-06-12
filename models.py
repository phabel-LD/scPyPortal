from sqlmodel import SQLModel, Field, Column, JSON, Relationship
from typing import Optional, List, Dict, Any
from datetime import datetime, timedelta
from pydantic import BaseModel
import json

class DatasetCacheBase(SQLModel):
    """Base model for dataset caching"""
    file_hash: str = Field(index=True, unique=True, max_length=64)
    original_filename: str
    file_size_mb: float
    upload_date: datetime = Field(default_factory=datetime.now)
    
    # Dataset metrics
    n_cells: int
    n_genes: int
    n_hvgs: Optional[int] = None
    pca_variance: List[float] = Field(sa_column=Column(JSON))
    metadata_summary: Dict[str, Any] = Field(sa_column=Column(JSON))
    
    # Processing status
    is_processed: bool = Field(default=False)
    processing_time_sec: Optional[float] = None
    
    # File paths
    processed_data_path: str  # Relative path to cached .h5ad

class DatasetCache(DatasetCacheBase, table=True):
    """Database table for dataset cache"""
    id: Optional[int] = Field(default=None, primary_key=True)
    
    # Relationships
    analysis_sessions: List["AnalysisSession"] = Relationship(back_populates="dataset")

class AnalysisSession(SQLModel, table=True):
    """Track user analysis sessions"""
    id: Optional[int] = Field(default=None, primary_key=True)
    session_id: str = Field(index=True)
    dataset_id: int = Field(foreign_key="datasetcache.id")
    start_time: datetime = Field(default_factory=datetime.now)
    end_time: Optional[datetime] = None
    actions_performed: List[str] = Field(sa_column=Column(JSON))
    
    # Relationships
    dataset: DatasetCache = Relationship(back_populates="analysis_sessions")

class PrecomputedResults(SQLModel, table=True):
    """Cache expensive computation results"""
    id: Optional[int] = Field(default=None, primary_key=True)
    dataset_hash: str = Field(index=True)
    computation_type: str  # e.g., "de_results", "paga", "dpt"
    parameters: Dict[str, Any] = Field(sa_column=Column(JSON))
    results: Dict[str, Any] = Field(sa_column=Column(JSON))
    computed_at: datetime = Field(default_factory=datetime.now)

# Pydantic models for API responses
class DatasetSummary(BaseModel):
    filename: str
    n_cells: int
    n_genes: int
    metadata_fields: List[str]
    is_cached: bool
    cache_age_days: float

class ProcessingParameters(BaseModel):
    n_neighbors: int = 15
    resolution: float = 0.5
    n_pcs: int = 50
    perplexity: int = 30
    random_state: int = 42