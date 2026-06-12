import sqlite3
import pickle
import zlib
import hashlib
from datetime import datetime, timedelta
from typing import Any, Optional, Dict
import threading
import asyncio
from contextlib import contextmanager

class IntelligentCache:
    """
    Advanced caching system with compression, expiration, and memory management
    """
    
    def __init__(self, db_path: str = "scpyportal_cache.db", max_cache_size_gb: float = 5.0):
        self.db_path = db_path
        self.max_cache_size = max_cache_size_gb * 1024 ** 3  # Convert to bytes
        self._init_database()
        self._cleanup_lock = threading.Lock()
    
    def _init_database(self):
        """Initialize cache database with optimized settings"""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS cache (
                    key TEXT PRIMARY KEY,
                    value BLOB,
                    compressed BOOLEAN DEFAULT 0,
                    size_bytes INTEGER,
                    access_count INTEGER DEFAULT 0,
                    last_accessed TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    expires_at TIMESTAMP,
                    computation_type TEXT,
                    parameters_hash TEXT
                )
            """)
            
            # Create indexes for performance
            conn.execute("CREATE INDEX IF NOT EXISTS idx_expires ON cache(expires_at)")
            conn.execute("CREATE INDEX IF NOT EXISTS idx_computation ON cache(computation_type)")
            conn.execute("CREATE INDEX IF NOT EXISTS idx_last_accessed ON cache(last_accessed)")
            
            # Enable WAL mode for better concurrent performance
            conn.execute("PRAGMA journal_mode=WAL")
            conn.execute("PRAGMA synchronous=NORMAL")
            conn.execute("PRAGMA cache_size=-64000")  # 64MB cache
    
    def _compute_key(self, computation_type: str, parameters: Dict) -> str:
        """Compute deterministic cache key"""
        param_str = str(sorted(parameters.items()))
        return hashlib.sha256(f"{computation_type}:{param_str}".encode()).hexdigest()
    
    def _compress_data(self, data: Any) -> bytes:
        """Compress data for storage"""
        return zlib.compress(pickle.dumps(data), level=6)
    
    def _decompress_data(self, compressed_data: bytes) -> Any:
        """Decompress stored data"""
        return pickle.loads(zlib.decompress(compressed_data))
    
    async def get(self, computation_type: str, parameters: Dict, 
                 max_age_hours: Optional[int] = None) -> Optional[Any]:
        """Retrieve item from cache with async support"""
        key = self._compute_key(computation_type, parameters)
        
        def _sync_get():
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.execute("""
                    SELECT value, compressed, expires_at 
                    FROM cache 
                    WHERE key = ? AND (expires_at IS NULL OR expires_at > ?)
                """, (key, datetime.now()))
                
                row = cursor.fetchone()
                if row:
                    # Update access statistics
                    conn.execute("""
                        UPDATE cache 
                        SET access_count = access_count + 1, last_accessed = ?
                        WHERE key = ?
                    """, (datetime.now(), key))
                    
                    value, compressed, expires_at = row
                    
                    if compressed:
                        return self._decompress_data(value)
                    else:
                        return pickle.loads(value)
                return None
        
        # Run in thread pool to avoid blocking
        loop = asyncio.get_event_loop()
        return await loop.run_in_executor(None, _sync_get)
    
    async def set(self, computation_type: str, parameters: Dict, value: Any,
                 expire_hours: Optional[int] = 24, compress: bool = True):
        """Store item in cache with async support"""
        key = self._compute_key(computation_type, parameters)
        
        def _sync_set():
            # Serialize value
            if compress:
                serialized = self._compress_data(value)
                size = len(serialized)
            else:
                serialized = pickle.dumps(value)
                size = len(serialized)
            
            # Calculate expiration
            expires_at = None
            if expire_hours:
                expires_at = datetime.now() + timedelta(hours=expire_hours)
            
            with sqlite3.connect(self.db_path) as conn:
                conn.execute("""
                    INSERT OR REPLACE INTO cache 
                    (key, value, compressed, size_bytes, computation_type, 
                     parameters_hash, expires_at, created_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, (key, serialized, compress, size, computation_type, 
                      str(hash(frozenset(parameters.items()))), expires_at, datetime.now()))
            
            # Check cache size and cleanup if needed
            self._enforce_size_limits()
        
        loop = asyncio.get_event_loop()
        await loop.run_in_executor(None, _sync_set)
    
    def _enforce_size_limits(self):
        """Enforce cache size limits by removing least recently used items"""
        with self._cleanup_lock:
            with sqlite3.connect(self.db_path) as conn:
                # Get current cache size
                total_size = conn.execute("SELECT SUM(size_bytes) FROM cache").fetchone()[0] or 0
                
                if total_size > self.max_cache_size:
                    # Remove least recently used items until under limit
                    items_to_remove = []
                    cursor = conn.execute("""
                        SELECT key, size_bytes 
                        FROM cache 
                        ORDER BY last_accessed ASC
                    """)
                    
                    current_size = total_size
                    for key, size in cursor:
                        if current_size <= self.max_cache_size * 0.8:  # Leave some headroom
                            break
                        
                        items_to_remove.append(key)
                        current_size -= size
                    
                    # Remove selected items
                    if items_to_remove:
                        placeholders = ','.join('?' * len(items_to_remove))
                        conn.execute(f"DELETE FROM cache WHERE key IN ({placeholders})", items_to_remove)
    
    async def cleanup_expired(self):
        """Remove expired cache entries"""
        def _sync_cleanup():
            with sqlite3.connect(self.db_path) as conn:
                conn.execute("DELETE FROM cache WHERE expires_at < ?", (datetime.now(),))
        
        loop = asyncio.get_event_loop()
        await loop.run_in_executor(None, _sync_cleanup)
    
    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics"""
        with sqlite3.connect(self.db_path) as conn:
            stats = {}
            
            # Basic stats
            stats['total_entries'] = conn.execute("SELECT COUNT(*) FROM cache").fetchone()[0]
            stats['total_size_mb'] = conn.execute("SELECT SUM(size_bytes) FROM cache").fetchone()[0] or 0
            stats['total_size_mb'] = stats['total_size_mb'] / (1024 ** 2)
            
            # Popular computations
            popular = conn.execute("""
                SELECT computation_type, COUNT(*), SUM(access_count)
                FROM cache 
                GROUP BY computation_type 
                ORDER BY SUM(access_count) DESC
            """).fetchall()
            
            stats['popular_computations'] = [
                {'type': row[0], 'count': row[1], 'total_accesses': row[2]} 
                for row in popular
            ]
            
            return stats