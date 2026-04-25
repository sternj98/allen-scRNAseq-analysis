"""
Data loading utilities for Allen Brain Map scRNA-seq data.

This module provides functions to load Allen Institute Cell Types Database
scRNA-seq data, including expression matrices and cell metadata.

Allen VISp data format:
- Expression matrix: genes (rows) × samples (columns)
- Metadata: sample_name column links to expression column headers

Recommended usage
-----------------
Use load_allen_filtered() instead of the two-step load_allen_data() +
create_anndata() workflow. It filters cells using the metadata file (9 MB)
before touching the expression CSV (1.5 GB), reducing peak memory ~10x.

    adata = load_allen_filtered(
        metadata_path=metadata_path,
        expression_path=expression_path,
        target_class="GABAergic",
        target_layers=["L2/3"],
    )
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Tuple, Optional, List
import scanpy as sc


def load_allen_expression(
    expression_path: str,
    cell_ids: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Load Allen Brain Map expression matrix.

    Allen format: rows = genes, columns = samples (cells)
    First row contains sample_name identifiers
    First column contains gene identifiers

    Parameters
    ----------
    expression_path : str
        Path to expression matrix CSV (e.g., mouse_VISp_2018-06-14_exon-matrix.csv)
    cell_ids : list of str, optional
        If provided, load only these cell columns. Dramatically reduces memory
        when loading a filtered subset. IDs not found in the CSV are silently
        skipped with a warning.

    Returns
    -------
    pd.DataFrame
        Expression matrix transposed to (cells × genes) format
        Index = sample_name, Columns = gene symbols
    """
    print(f"Loading expression matrix from {expression_path}...")

    if cell_ids is not None:
        # Read only the header to learn column names, then select target columns
        all_cols = pd.read_csv(expression_path, nrows=0).columns.tolist()
        gene_col = all_cols[0]          # 'Unnamed: 0' — the row index column
        file_cell_set = set(all_cols[1:])
        cols_to_load = [gene_col] + [c for c in cell_ids if c in file_cell_set]

        missing = len(cell_ids) - (len(cols_to_load) - 1)
        if missing > 0:
            print(f"  Warning: {missing} requested cell IDs not found in expression file")

        print(f"  Reading {len(cols_to_load) - 1} of {len(all_cols) - 1} cell columns...")
        expression = pd.read_csv(expression_path, usecols=cols_to_load, index_col=0)
    else:
        expression = pd.read_csv(expression_path, index_col=0)

    # Expression is genes × samples, transpose to samples × genes
    expression = expression.T

    print(f"Loaded expression: {expression.shape[0]} cells × {expression.shape[1]} genes")
    return expression


def load_allen_metadata(
    metadata_path: str,
    sample_col: str = "sample_name"
) -> pd.DataFrame:
    """
    Load Allen Brain Map cell metadata.

    Parameters
    ----------
    metadata_path : str
        Path to metadata CSV (e.g., mouse_VISp_2018-06-14_samples-columns.csv)
    sample_col : str
        Column name containing sample identifiers (default: 'sample_name')

    Returns
    -------
    pd.DataFrame
        Cell metadata indexed by sample_name
        Key columns: class, subclass, cluster, brain_subregion, etc.
    """
    print(f"Loading metadata from {metadata_path}...")

    metadata = pd.read_csv(metadata_path)

    # Set sample_name as index
    if sample_col in metadata.columns:
        metadata = metadata.set_index(sample_col)
    else:
        raise ValueError(f"Column '{sample_col}' not found. Available: {metadata.columns.tolist()}")

    print(f"Loaded metadata: {metadata.shape[0]} cells × {metadata.shape[1]} columns")
    print(f"Key columns: {[c for c in ['class', 'subclass', 'cluster'] if c in metadata.columns]}")

    return metadata


def load_allen_data(
    metadata_path: str,
    expression_path: str,
    sample_col: str = "sample_name"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load Allen Brain Map scRNA-seq data from CSV files.

    Parameters
    ----------
    metadata_path : str
        Path to samples-columns.csv metadata file
    expression_path : str
        Path to exon-matrix.csv expression file
    sample_col : str
        Column name for sample identifiers in metadata

    Returns
    -------
    metadata : pd.DataFrame
        Cell metadata indexed by sample_name
    expression : pd.DataFrame
        Expression matrix (cells × genes) indexed by sample_name

    Example
    -------
    >>> metadata, expression = load_allen_data(
    ...     'mouse_VISp_2018-06-14_samples-columns.csv',
    ...     'mouse_VISp_2018-06-14_exon-matrix.csv'
    ... )
    """
    metadata = load_allen_metadata(metadata_path, sample_col=sample_col)
    expression = load_allen_expression(expression_path)

    # Verify alignment
    common_samples = metadata.index.intersection(expression.index)
    n_meta_only = len(metadata.index.difference(expression.index))
    n_expr_only = len(expression.index.difference(metadata.index))

    if n_meta_only > 0 or n_expr_only > 0:
        print(f"Note: {n_meta_only} samples in metadata only, {n_expr_only} in expression only")

    print(f"Common samples: {len(common_samples)}")

    return metadata, expression


def create_anndata(
    metadata: pd.DataFrame,
    expression: pd.DataFrame,
    layer_name: Optional[str] = None,
    genes_path: Optional[str] = None
) -> sc.AnnData:
    """
    Create AnnData object from metadata and expression DataFrames.

    Parameters
    ----------
    metadata : pd.DataFrame
        Cell metadata (indexed by sample_name)
    expression : pd.DataFrame
        Expression matrix (cells × genes, indexed by sample_name)
    layer_name : str, optional
        Name for expression layer (e.g., 'counts', 'tpm')
    genes_path : str, optional
        Path to genes-rows.csv for gene annotations

    Returns
    -------
    sc.AnnData
        AnnData object with:
        - .X: expression matrix
        - .obs: cell metadata
        - .var: gene annotations (if genes_path provided)
        - .obs_names: sample_name identifiers
        - .var_names: gene symbols
    """
    # Align cells between metadata and expression
    # Convert indices to strings first for proper matching
    metadata = metadata.copy()
    metadata.index = metadata.index.astype(str)
    expression = expression.copy()
    expression.index = expression.index.astype(str)
    expression.columns = expression.columns.astype(str)

    common_cells = metadata.index.intersection(expression.index)
    print(f"Creating AnnData with {len(common_cells)} cells")

    metadata_aligned = metadata.loc[common_cells].copy()
    expression_aligned = expression.loc[common_cells].copy()

    # Ensure obs index is a proper string Index
    metadata_aligned.index = pd.Index(metadata_aligned.index.astype(str))

    # Create AnnData object with explicit string indices
    adata = sc.AnnData(
        X=expression_aligned.values.astype(np.float32),
        obs=metadata_aligned,
    )
    adata.var_names = pd.Index(expression_aligned.columns.astype(str))

    # Add gene annotations if provided
    # Allen genes-rows.csv links to the expression matrix via gene_entrez_id,
    # which matches the numeric row index of the expression CSV.
    if genes_path:
        genes_df = pd.read_csv(genes_path)
        if 'gene_entrez_id' in genes_df.columns:
            genes_df = genes_df.set_index('gene_entrez_id')
        elif 'gene' in genes_df.columns:
            genes_df = genes_df.set_index('gene')
        else:
            genes_df = genes_df.set_index(genes_df.columns[0])
        genes_df.index = genes_df.index.astype(str)
        adata.var = genes_df.reindex(adata.var_names)

    if layer_name:
        adata.layers[layer_name] = adata.X.copy()

    # Print summary of cell types
    if 'class' in adata.obs.columns:
        print(f"\nCell classes: {adata.obs['class'].nunique()}")
        print(adata.obs['class'].value_counts())

    return adata


def _load_obsm_from_paths(
    obsm_paths: Dict[str, str],
    cell_ids: List[str],
) -> Dict[str, np.ndarray]:
    """
    Load pre-computed cell embeddings from files and align to cell_ids order.

    Each file must be one of:
      - CSV  (.csv)  — index column = sample_name, remaining columns = dimensions
      - NumPy (.npy) — shape (n_cells, n_dims), rows already aligned to cell_ids

    Returns a dict suitable for direct assignment to adata.obsm.
    Cells missing from a CSV file receive NaN rows; .npy arrays must already
    match len(cell_ids) exactly.
    """
    result: Dict[str, np.ndarray] = {}
    cell_index = pd.Index(cell_ids)

    for key, path in obsm_paths.items():
        p = Path(path)
        if not p.exists():
            print(f"  obsm[{key!r}]: file not found — {path}")
            continue

        if p.suffix.lower() == ".npy":
            arr = np.load(path).astype(np.float32)
            if arr.shape[0] != len(cell_ids):
                raise ValueError(
                    f"obsm[{key!r}]: .npy array has {arr.shape[0]} rows "
                    f"but {len(cell_ids)} cells were selected"
                )
            result[key] = arr

        else:
            df = pd.read_csv(path, index_col=0)
            df.index = df.index.astype(str)
            missing = cell_index.difference(df.index)
            if len(missing):
                print(f"  obsm[{key!r}]: {len(missing)} cells not in file — filled with NaN")
            result[key] = df.reindex(cell_index).values.astype(np.float32)

        print(f"  obsm[{key!r}]: loaded shape {result[key].shape}")

    return result


def load_allen_filtered(
    metadata_path: str,
    expression_path: str,
    sample_col: str = "sample_name",
    target_class: Optional[str] = "GABAergic",
    class_col: str = "class",
    target_layers: Optional[List[str]] = None,
    layer_col: str = "brain_subregion",
    genes_path: Optional[str] = None,
    layer_name: Optional[str] = None,
    obsm_paths: Optional[Dict[str, str]] = None,
) -> sc.AnnData:
    """
    Load Allen data filtered to a cell subset without reading the full expression matrix.

    Loads the metadata CSV first (fast, ~9 MB), applies cell-level filters, then
    reads only the matching columns from the expression CSV. Compared to loading
    the full matrix and filtering afterwards, peak memory drops ~10x for typical
    Layer 2/3 subsets (1173 of 15413 cells).

    Parameters
    ----------
    metadata_path : str
        Path to samples-columns.csv
    expression_path : str
        Path to exon-matrix.csv
    sample_col : str
        Column in metadata that matches expression CSV column headers
    target_class : str, optional
        Keep only cells with this class label (e.g., 'GABAergic').
        Pass None to skip class filtering.
    class_col : str
        Metadata column for cell class
    target_layers : list of str, optional
        Keep only cells whose layer_col value is in this list
        (e.g., ['L2/3']). Pass None to skip layer filtering.
        Allen brain_subregion values are 'L2/3', 'L5', 'L6', etc.
    layer_col : str
        Metadata column for cortical layer
    genes_path : str, optional
        Path to genes-rows.csv for gene annotations
    layer_name : str, optional
        Name to store expression in adata.layers (e.g., 'counts')
    obsm_paths : dict of str -> str, optional
        Pre-computed cell embeddings to load into adata.obsm. Keys become
        the obsm slot names (e.g., 'X_umap', 'X_tsne', 'X_pca'); values
        are file paths. Supported formats:
          - CSV  (.csv): index = sample_name, columns = embedding dimensions
          - NumPy (.npy): array of shape (n_selected_cells, n_dims)

        Example::

            obsm_paths={
                'X_tsne': '/path/to/tsne_coords.csv',
                'X_umap': '/path/to/umap_coords.npy',
            }

    Returns
    -------
    sc.AnnData
        AnnData filtered to the requested cell subset, with obsm populated
        if obsm_paths was provided.

    Example
    -------
    >>> adata = load_allen_filtered(
    ...     metadata_path=metadata_path,
    ...     expression_path=expression_path,
    ...     target_class="GABAergic",
    ...     target_layers=["L2/3"],
    ... )
    """
    if target_layers is None:
        target_layers = ["L2/3"]

    # ── Step 1: load metadata (fast) ──────────────────────────────────────────
    metadata = load_allen_metadata(metadata_path, sample_col=sample_col)

    # ── Step 2: apply cell-level filters on metadata ──────────────────────────
    mask = pd.Series(True, index=metadata.index)

    if target_class is not None and class_col in metadata.columns:
        mask &= metadata[class_col] == target_class

    if target_layers and layer_col in metadata.columns:
        mask &= metadata[layer_col].isin(target_layers)

    target_ids = metadata.index[mask].tolist()
    n_total = len(metadata)
    print(
        f"Metadata filter: {len(target_ids)} of {n_total} cells selected "
        f"(class={target_class!r}, layers={target_layers})"
    )

    if len(target_ids) == 0:
        avail_class = metadata[class_col].unique().tolist() if class_col in metadata.columns else "N/A"
        avail_layer = metadata[layer_col].unique().tolist() if layer_col in metadata.columns else "N/A"
        raise ValueError(
            f"No cells matched the filter criteria.\n"
            f"  Available {class_col}: {avail_class}\n"
            f"  Available {layer_col}: {avail_layer}"
        )

    # ── Step 3: load only the target cell columns from expression CSV ─────────
    expression = load_allen_expression(expression_path, cell_ids=target_ids)

    # ── Step 4: build AnnData ─────────────────────────────────────────────────
    metadata_filtered = metadata.loc[metadata.index.intersection(expression.index)]
    adata = create_anndata(metadata_filtered, expression, layer_name=layer_name, genes_path=genes_path)

    # ── Step 5: populate obsm from external embedding files ───────────────────
    if obsm_paths:
        print("Loading obsm embeddings...")
        # Use the final cell order from adata, which may differ from target_ids
        # after create_anndata aligns metadata and expression
        embeddings = _load_obsm_from_paths(obsm_paths, adata.obs_names.tolist())
        for key, arr in embeddings.items():
            adata.obsm[key] = arr

    return adata


def filter_by_quality(
    adata: sc.AnnData,
    min_genes: int = 500,
    min_counts: int = 1000,
    max_mito_pct: float = 20.0,
    mito_prefix: str = "mt-"
) -> sc.AnnData:
    """
    Apply basic quality control filters to scRNA-seq data.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    min_genes : int
        Minimum number of genes per cell
    min_counts : int
        Minimum total counts per cell
    max_mito_pct : float
        Maximum mitochondrial gene percentage
    mito_prefix : str
        Prefix for mitochondrial genes

    Returns
    -------
    sc.AnnData
        Filtered AnnData object
    """
    n_cells_before = adata.n_obs

    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith(mito_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Apply filters
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, min_counts=min_counts)

    if 'pct_counts_mt' in adata.obs.columns:
        adata = adata[adata.obs['pct_counts_mt'] < max_mito_pct, :]

    n_cells_after = adata.n_obs
    print(f"Filtered {n_cells_before - n_cells_after} cells ({n_cells_after} remaining)")

    return adata
