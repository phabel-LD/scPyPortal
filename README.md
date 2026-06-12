scPyPortal v1 – Interactive Single‑Cell Analysis Platform
================================================================

scPyPortal is a full‑featured, interactive web application for single‑cell RNA‑seq (scRNA‑seq) and multi‑omics data analysis, built with Shiny for Python, Scanpy, Muon, and Plotly. It provides a point‑and‑click interface to perform quality control, clustering, differential expression, trajectory inference, and multi‑omics integration – all without writing a single line of code.

Table of Contents
-----------------
1. Key Features
2. Architecture Overview
3. Installation
4. Quick Start
5. Detailed User Guide
6. Module Reference
7. Configuration & Caching
8. Troubleshooting
9. Extending the Application
10. Citation & Acknowledgments

1. Key Features
---------------
Feature                      | Description
-----------------------------|---------------------------------------------------------------
File Upload                  | Supports .h5ad (AnnData) and .h5mu (MuData multi‑omics).
Quality Control              | Interactive histograms for genes/cell, UMI counts, mitochondrial percentage.
Clustering                   | PCA, kNN graph, UMAP, Leiden clustering with adjustable resolution.
Differential Expression      | T‑test, Wilcoxon, logistic regression; volcano plots, marker tables, gene expression violin plots.
Trajectory Inference         | Diffusion maps and DPT pseudotime; root cell selection; gene expression trends along pseudotime.
Multi‑omics Integration      | WNN (Weighted Nearest Neighbors) and NMF‑based integration of scRNA+scATAC data.
Smart Caching                | SQLite‑based cache with compression, expiration, and LRU eviction.
Export                       | Download the fully processed AnnData object (.h5ad) with all embeddings, clusters, and results.
Asynchronous Processing      | Heavy computations run in background – UI never freezes.

2. Architecture Overview
------------------------
The application follows a modular reactive design using Shiny for Python.

Directory structure:
scPyPortal/
├── app.py                    # Main UI + server, global reactive state
├── config.py                 # Paths, default parameters, directory setup
├── models.py                 # SQLModel database schemas (cache, sessions)
├── modules/                  # Independent analysis modules
│   ├── preprocessing.py      # QC metrics and filtering
│   ├── clustering.py         # PCA, UMAP, Leiden, silhouette
│   ├── differential_expression.py
│   ├── trajectory.py
│   ├── multiomics.py
│   └── export.py
├── utils/                    # Helper libraries
│   ├── cache_manager.py      # High‑level caching with SQLModel
│   ├── advanced_cache.py     # Compressed, TTL, LRU cache
│   ├── data_loader.py        # Load from cache or disk
│   ├── optimized_loader.py   # Backed / chunked loading for huge datasets
│   ├── validators.py         # File and parameter validation
│   ├── visualization.py      # Plotly helper functions
│   └── hashing.py            # SHA256 file hashing
└── cache/                    # Auto‑created directories for cached data
    ├── datasets/             # Stored .h5ad files
    └── temp/

Data Flow (simplified):
1. Upload → app.py computes hash, checks DatasetCache, stores file in cache/datasets/.
2. User interaction → Module buttons trigger reactive.Effect.
3. Processing → Modules copy current_adata(), run Scanpy/Muon methods.
4. Caching → Expensive results are stored via IntelligentCache (compressed, keyed by parameters).
5. Visualization → Matplotlib/Plotly figures are rendered reactively.
6. Export → Current AnnData written to a temporary H5AD and streamed to browser.

3. Installation
---------------
Prerequisites:
- Python 3.10 or 3.11
- At least 8 GB RAM (16+ recommended for large datasets)
- Conda (optional but recommended)

Step‑by‑step:

1. Clone the repository
   git clone https://github.com/yourusername/scPyPortal.git
   cd scPyPortal

2. Create and activate a virtual environment
   # Using conda
   conda create -n scpyportal python=3.10 -y
   conda activate scpyportal

   # Or using venv
   python -m venv scpyportal_env
   source scpyportal_env/bin/activate   # Linux/Mac
   scpyportal_env\Scripts\activate      # Windows

3. Install dependencies
   pip install -r requirements.txt

   Note: Some packages (e.g., muon) may require additional system libraries.
   On Ubuntu: sudo apt-get install libhdf5-dev

4. Run the app
   shiny run app.py

   Open your browser at http://127.0.0.1:8000.

Optional: Set environment variables for production (e.g., SHINY_APP_HOST=0.0.0.0 to serve on network).

4. Quick Start
--------------
1. Upload a dataset – Click "Upload Single‑Cell Data" and choose an .h5ad file (e.g., the classic pbmc3k.h5ad from Scanpy).
2. Go to Quality Control – Click "Run QC" to see histograms of genes/cell, UMI counts, and mitochondrial percentage.
3. Go to Clustering – Adjust Resolution (try 0.5) and click "Run Clustering". Explore UMAP and cluster sizes.
4. Go to Differential Expression – Select a grouping variable (e.g., leiden), choose two clusters, and run DE. View volcano plot and table.
5. Go to Trajectory Analysis – Choose a root cluster, run DPT, and visualise pseudotime.
6. Export – Download your fully processed AnnData file for further analysis in R/Seurat or other tools.

5. Detailed User Guide
----------------------

Sidebar (Global Controls)
  * File upload: Load .h5ad or .h5mu. The dataset is cached automatically.
  * Dataset Info: Shows filename, number of cells/genes, metadata fields.
  * QC Status: Displays mean genes/cell and counts after QC.
  * QC Plots: Histograms of QC metrics (only shown after QC is run).
  * Analysis Controls: "Run Preprocessing" (basic QC), "Clear Cache" (removes cached files).

Module Tabs

1. Quality Control
   Parameters: Min/max genes per cell, max mitochondrial % (filtering not yet implemented – only visual QC).
   Run QC – calculates metrics and displays histograms.
   Output: Genes/cell, UMI counts, mito% distributions.

2. Clustering
   Parameters: Number of PCs, kNN neighbours, Leiden resolution, distance metric.
   Run Clustering – performs PCA, neighbours, UMAP, Leiden.
   Results:
     - Number of clusters + Silhouette score (based on PCA).
     - UMAP plot (colour by cluster or any categorical variable).
     - Bar plot of cluster sizes.

3. Differential Expression
   Group By – choose a categorical column from adata.obs.
   Group 1 / Group 2 – select two groups to compare.
   Test method – t‑test (default), Wilcoxon, logistic regression.
   Filters: min log2 fold change, max adjusted p‑value, min cells per group.
   Run DE Analysis – computes markers using scanpy.tl.rank_genes_groups.
   Find All Markers – computes markers for every group against the rest (stored in memory).
   Tabs:
     - Volcano Plot – log2FC vs -log10(adj. p‑value).
     - DE Results – interactive datatable with gene, log2FC, p‑value, score.
     - Gene Expression – select a gene to see expression violin plot across groups.
     - All Markers – choose a cell type to see its top markers.

4. Trajectory Analysis
   Method – Diffusion Map (only dimensionality reduction) or DPT (pseudotime).
   Root Cell Grouping – pick a categorical column to define the root.
   Root Cluster/Type – the specific group that will be used as trajectory origin.
   Run Trajectory Analysis:
     - For DPT: sets iroot to a random cell from the chosen group, computes diffusion map, then DPT.
   Output tabs:
     - Diffusion Map – scatter plot of DC1 vs DC2.
     - Pseudotime – overlay of pseudotime on UMAP (or diffusion map if UMAP absent).
     - Gene Trends – select a gene to see its expression smoothed along pseudotime.

5. Multi‑omics Integration (requires .h5mu file)
   Method – WNN (Weighted Nearest Neighbors) or NMF.
   Modalities to Integrate – select which assays (e.g., rna, atac) to combine.
   Run Integration – computes joint embedding.
   Output tabs:
     - Integration UMAP – combined UMAP coloured by chosen metadata.
     - Modality Weights (WNN only) – bar plot of modality contributions.
     - Cross‑modality – scatter plot of gene expression vs peak accessibility, with correlation.
   Note: MOFA and DIABLO are placeholders for future versions.

6. Export Results
   Export filename – specify name for the downloaded file.
   Download .h5ad File – writes the current AnnData object (including PCA, UMAP, clusters, DE results, pseudotime) to a temporary file and sends it to the browser.

6. Module Reference
-------------------
Each module (preprocessing, clustering, differential_expression, trajectory, multiomics, export) exposes two functions:

- create_<module>_ui(id: str) -> ui.Component
  Returns the Shiny UI elements for that module, using ui.ns(id) to isolate input/output IDs.

- create_<module>_server(id, current_adata, is_processing) -> dict
  Defines reactive logic, returns a dictionary with internal reactive values (e.g., {'adata_processed': ..., 'clustering_results': ...}).
  Parameters:
    - id : module identifier (must match UI call)
    - current_adata : reactive.Value containing the AnnData/MuData object
    - is_processing : reactive.Value to disable buttons during computation

7. Configuration & Caching
--------------------------
All settings are in config.py:
  CACHE_DIR = BASE_DIR / "cache"
  DATABASE_URL = "sqlite:///./scpyportal.db"
  DEFAULT_N_NEIGHBORS = 15
  DEFAULT_RESOLUTION = 0.5
  MAX_FILE_SIZE_MB = 500
  CACHE_EXPIRY_DAYS = 30

Two caching layers:
  * DatasetCache (SQLModel) – stores entire .h5ad files and metadata; avoids re‑uploading.
  * IntelligentCache (advanced_cache.py) – stores serialized results of computations (DE, clustering, etc.) with:
      - Key = hash(computation_type + sorted parameters)
      - Value = compressed (zlib) pickle
      - Expiration (TTL) and LRU eviction when total cache exceeds max_cache_size_gb (default 5 GB).
      - Async get/set using asyncio + run_in_executor.

8. Troubleshooting
------------------
Problem: App does not start – "ModuleNotFoundError"
  Solution: Ensure all dependencies are installed: pip install -r requirements.txt

Problem: File upload fails with "Cannot read .h5ad"
  Solution: Verify the file is a valid AnnData HDF5 file (use scanpy.read_h5ad from command line).

Problem: Memory error with large dataset (>2GB)
  Solution: utils/optimized_loader.py automatically switches to backed mode. Ensure you have enough disk space for cache.

Problem: Clustering button does nothing
  Check browser console (F12) for errors. Ensure that QC or preprocessing has been run at least once (to have log-normalised data).

Problem: Multi‑omics tab shows nothing
  You must upload a .h5mu (MuData) file, not a regular .h5ad.

Problem: Downloaded .h5ad file is corrupted
  This can happen on very large files; try setting a smaller export name without special characters. The app writes to a temporary file then streams – check disk space.

9. Extending the Application
----------------------------
To add a new analysis module:
1. Create a new file modules/my_module.py
2. Implement create_my_module_ui(id) and create_my_module_server(id, current_adata, is_processing).
3. Import the module in app.py.
4. Add a new ui.nav_panel to create_main_interface().
5. Instantiate the server inside server() function, passing the reactive values.

To add a new caching backend:
  Subclass IntelligentCache and override _compress_data/_decompress_data, or add new storage (e.g., Redis) by implementing similar async methods.

To support additional file formats (e.g., .loom, .zarr):
  Extend data_loader.py and optimized_loader.py with new load methods and update the file upload accept list.

10. Citation & Acknowledgments
------------------------------
scPyPortal v1 was developed as a final project for a Data Analysis diploma course.
Author: Phabel Antonio Lopez Delgado
Contact: phabel2001@gmail.com | phabel@lcg.unam.mx

This software uses the following open‑source libraries:
  - Shiny for Python (Posit)
  - Scanpy (Theis Lab)
  - Muon (scverse)
  - Plotly
  - SQLModel (FastAPI team)

License: MIT