"""
Cell filtering utilities for subsetting to specific cell populations.

This module provides functions to filter cells based on Allen Institute taxonomy
annotations, focusing on Layer 2/3 inhibitory interneurons.

Allen VISp data columns:
- class: GABAergic, Glutamatergic, Non-Neuronal
- subclass: Pvalb, Sst, Vip, Lamp5, etc.
- cluster: Fine-grained cell type clusters
- brain_subregion: Contains layer info (e.g., VISp2/3, VISp5, etc.)

Caching
-------
All filter functions support disk caching of filtered AnnData objects.
Set a default cache directory once via set_cache_dir(), or pass cache_dir=
to any individual call. Each cache entry saves a .h5ad file and a .json
sidecar with filter metadata (function name, params, cell counts, timestamp).

    import src.cell_filter as cf
    cf.set_cache_dir(DATA_DIR / "cache")

    # subsequent calls load from cache automatically
    adata_inh = cf.filter_inhibitory(adata)
"""

import hashlib
import inspect
import json
import numpy as np
import pandas as pd
import scanpy as sc
from datetime import datetime
from functools import wraps
from pathlib import Path
from typing import List, Optional, Union


# Module-level default cache directory. Set via set_cache_dir() before calling
# filter functions, or pass cache_dir= to individual calls.
DEFAULT_CACHE_DIR: Optional[Path] = None


def set_cache_dir(path: Union[str, Path, None]) -> None:
    """Set the default cache directory used by all cached filter functions."""
    global DEFAULT_CACHE_DIR
    DEFAULT_CACHE_DIR = Path(path) if path is not None else None


def _obs_hash(adata: sc.AnnData) -> str:
    names = "".join(adata.obs_names.tolist())
    return hashlib.md5(names.encode()).hexdigest()


def _params_hash(params: dict) -> str:
    s = json.dumps(params, sort_keys=True, default=str)
    return hashlib.md5(s.encode()).hexdigest()[:8]


def cached_filter(fn):
    """
    Decorator that adds transparent disk caching to AnnData filter functions.

    Decorated functions gain two extra keyword arguments:
      cache_dir : str or Path, optional
          Directory for cache files. Overrides DEFAULT_CACHE_DIR for this call.
      use_cache : bool, default True
          Set to False to bypass cache entirely (re-runs filter and overwrites).

    Cache files are named <fn>_<input_hash>_<params_hash>.h5ad with a
    companion .json metadata sidecar written alongside.
    """
    sig = inspect.signature(fn)

    @wraps(fn)
    def wrapper(adata, *args, cache_dir=None, use_cache=True, **kwargs):
        effective_cache_dir = Path(cache_dir) if cache_dir is not None else DEFAULT_CACHE_DIR

        if not use_cache or effective_cache_dir is None:
            return fn(adata, *args, **kwargs)

        # Resolve full param dict (excluding adata)
        bound = sig.bind(adata, *args, **kwargs)
        bound.apply_defaults()
        params = {k: v for k, v in bound.arguments.items() if k != "adata"}

        obs_h = _obs_hash(adata)
        param_h = _params_hash(params)
        cache_key = f"{fn.__name__}_{obs_h[:8]}_{param_h}"

        cache_path = effective_cache_dir / f"{cache_key}.h5ad"
        meta_path = effective_cache_dir / f"{cache_key}.json"

        if cache_path.exists():
            print(f"[cache] Hit — loading {fn.__name__} from {cache_path.name}")
            return sc.read_h5ad(cache_path)

        result = fn(adata, *args, **kwargs)

        effective_cache_dir.mkdir(parents=True, exist_ok=True)
        metadata = {
            "function": fn.__name__,
            "timestamp": datetime.now().isoformat(),
            "cache_key": cache_key,
            "n_input_cells": adata.n_obs,
            "n_input_genes": adata.n_vars,
            "n_output_cells": result.n_obs,
            "n_output_genes": result.n_vars,
            "params": {k: str(v) for k, v in params.items()},
        }
        result.uns["_cache_metadata"] = metadata
        result.write_h5ad(cache_path)
        with open(meta_path, "w") as f:
            json.dump(metadata, f, indent=2)
        print(f"[cache] Saved {fn.__name__} result → {cache_path.name}")

        return result

    return wrapper


@cached_filter
def filter_inhibitory(
    adata: sc.AnnData,
    class_col: str = "class",
    inhibitory_class: str = "GABAergic"
) -> sc.AnnData:
    """
    Filter cells to inhibitory (GABAergic) neurons.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    class_col : str
        Column name for cell class (default: 'class')
    inhibitory_class : str
        Class label for inhibitory neurons (default: 'GABAergic')

    Returns
    -------
    sc.AnnData
        Filtered AnnData containing only inhibitory neurons
    """
    mask = (adata.obs[class_col] == inhibitory_class).values
    n_cells = mask.sum()

    print(f"Filtering to {inhibitory_class} neurons: {n_cells} cells")

    if n_cells == 0:
        print(f"WARNING: No cells found. Available classes: {adata.obs[class_col].unique().tolist()}")

    return adata[mask, :].copy()


@cached_filter
def filter_by_layer(
    adata: sc.AnnData,
    target_layers: Union[str, List[str]],
    layer_col: str = "brain_subregion"
) -> sc.AnnData:
    """
    Filter cells by cortical layer using brain_subregion column.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    target_layers : str or List[str]
        Target layer patterns to match (e.g., '2/3', ['2/3', '4'])
        Will match substrings in brain_subregion column
    layer_col : str
        Column containing layer info (default: 'brain_subregion')

    Returns
    -------
    sc.AnnData
        Filtered AnnData

    Notes
    -----
    Allen brain_subregion format: 'VISp2/3', 'VISp5', 'VISp6a', etc.
    This function matches the layer portion (e.g., '2/3' matches 'VISp2/3')
    """
    if isinstance(target_layers, str):
        target_layers = [target_layers]

    # Build regex pattern to match any of the target layers
    # Match layer at end of string (e.g., VISp2/3 -> matches '2/3')
    pattern = '|'.join([f'{layer}$' for layer in target_layers])

    mask = adata.obs[layer_col].str.contains(pattern, regex=True, na=False).values
    n_cells = mask.sum()

    print(f"Filtering to layers {target_layers}: {n_cells} cells")

    if n_cells == 0:
        print(f"WARNING: No cells found. Available values: {adata.obs[layer_col].unique().tolist()}")

    return adata[mask, :].copy()


@cached_filter
def filter_layer_23_inhibitory(
    adata: sc.AnnData,
    class_col: str = "class",
    layer_col: str = "brain_subregion",
    inhibitory_class: str = "GABAergic",
    target_layers: Optional[List[str]] = None
) -> sc.AnnData:
    """
    Filter cells to Layer 2/3 inhibitory interneurons.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object with cell type annotations
    class_col : str
        Column for cell class (default: 'class')
    layer_col : str
        Column for layer annotation (default: 'brain_subregion')
    inhibitory_class : str
        Name of inhibitory class (default: 'GABAergic')
    target_layers : List[str], optional
        Layer patterns to match (default: ['2/3'])

    Returns
    -------
    sc.AnnData
        Filtered AnnData containing only Layer 2/3 inhibitory neurons

    Notes
    -----
    Allen taxonomy: class > subclass > cluster
    Allen brain_subregion format: 'VISp2/3', 'VISp5', etc.
    """
    if target_layers is None:
        target_layers = ['2/3']

    print(f"Filtering to {inhibitory_class} neurons in layers: {target_layers}")

    # Filter by cell class
    mask_class = (adata.obs[class_col] == inhibitory_class).values

    # Filter by layer (pattern match on brain_subregion)
    pattern = '|'.join([f'{layer}$' for layer in target_layers])
    mask_layer = adata.obs[layer_col].str.contains(pattern, regex=True, na=False).values

    # Combine filters
    mask_combined = mask_class & mask_layer

    n_cells = mask_combined.sum()
    print(f"Found {n_cells} Layer 2/3 inhibitory neurons")

    if n_cells == 0:
        print("WARNING: No cells matched filters.")
        print(f"  {class_col} values: {adata.obs[class_col].unique().tolist()}")
        print(f"  {layer_col} values: {adata.obs[layer_col].unique().tolist()}")

    return adata[mask_combined, :].copy()


@cached_filter
def filter_by_subclass(
    adata: sc.AnnData,
    target_subclasses: Optional[List[str]] = None,
    subclass_col: str = "subclass"
) -> sc.AnnData:
    """
    Filter cells to specific inhibitory subclasses.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    target_subclasses : List[str], optional
        Target subclass names (e.g., ['Pvalb', 'Sst', 'Vip', 'Lamp5'])
        If None, returns all cells unchanged
    subclass_col : str
        Column name for subclass annotation

    Returns
    -------
    sc.AnnData
        Filtered AnnData

    Notes
    -----
    Common Allen inhibitory subclasses:
    - Pvalb: Parvalbumin-expressing basket/chandelier cells
    - Sst: Somatostatin-expressing Martinotti cells
    - Vip: Vasoactive intestinal peptide-expressing
    - Lamp5: Neurogliaform and single bouquet cells
    - Sncg: Sparse interneurons
    - Serpinf1: Layer 1 interneurons
    """
    if target_subclasses is None:
        print("No subclass filter applied")
        return adata

    print(f"Filtering to subclasses: {target_subclasses}")

    mask = adata.obs[subclass_col].isin(target_subclasses).values
    n_cells = mask.sum()

    print(f"Found {n_cells} cells in target subclasses")

    if n_cells == 0:
        print(f"Available subclasses: {adata.obs[subclass_col].unique().tolist()}")

    return adata[mask, :].copy()


@cached_filter
def filter_by_cluster(
    adata: sc.AnnData,
    cluster_col: str = "cluster",
    min_cells_per_cluster: int = 10,
    target_clusters: Optional[List[str]] = None
) -> sc.AnnData:
    """
    Filter by cluster, optionally removing small clusters.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    cluster_col : str
        Column name for cluster annotation
    min_cells_per_cluster : int
        Minimum cells required per cluster (ignored if target_clusters specified)
    target_clusters : List[str], optional
        Specific clusters to keep

    Returns
    -------
    sc.AnnData
        Filtered AnnData
    """
    if target_clusters is not None:
        mask = adata.obs[cluster_col].isin(target_clusters).values
        print(f"Filtering to {len(target_clusters)} specified clusters: {mask.sum()} cells")
        return adata[mask, :].copy()

    # Filter by minimum cell count
    cluster_counts = adata.obs[cluster_col].value_counts()
    valid_clusters = cluster_counts[cluster_counts >= min_cells_per_cluster].index

    n_removed = len(cluster_counts) - len(valid_clusters)
    print(f"Removing {n_removed} clusters with < {min_cells_per_cluster} cells")

    mask = adata.obs[cluster_col].isin(valid_clusters).values
    return adata[mask, :].copy()


@cached_filter
def filter_high_quality(
    adata: sc.AnnData,
    exclude_no_class: bool = True,
    core_only: bool = False
) -> sc.AnnData:
    """
    Filter to high-quality cells based on Allen QC annotations.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    exclude_no_class : bool
        Remove cells with class='No Class' (low-quality/doublets)
    core_only : bool
        Keep only 'Core' cells (high-confidence assignments)

    Returns
    -------
    sc.AnnData
        Filtered AnnData

    Notes
    -----
    Allen provides:
    - class='No Class': Low-quality, doublets, or outliers
    - core_intermediate_call: 'Core' (>90/100 bootstrap trials) or 'Intermediate'
    """
    mask = np.ones(adata.n_obs, dtype=bool)

    if exclude_no_class and 'class' in adata.obs.columns:
        class_mask = (adata.obs['class'] != 'No Class').values
        n_removed = (~class_mask).sum()
        mask = mask & class_mask
        print(f"Excluding 'No Class' cells: {n_removed} removed")

    if core_only and 'core_intermediate_call' in adata.obs.columns:
        core_mask = (adata.obs['core_intermediate_call'] == 'Core').values
        mask = mask & core_mask
        print(f"Keeping only 'Core' cells: {mask.sum()} remaining")

    return adata[mask, :].copy()


def get_cell_type_summary(
    adata: sc.AnnData,
    groupby: str = "subclass"
) -> pd.DataFrame:
    """
    Generate summary statistics for cell type groups.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    groupby : str
        Column name to group by

    Returns
    -------
    pd.DataFrame
        Summary table with cell counts and percentages
    """
    counts = adata.obs[groupby].value_counts()
    percentages = (counts / len(adata) * 100).round(2)

    summary = pd.DataFrame({
        'n_cells': counts,
        'percentage': percentages
    })

    return summary.sort_values('n_cells', ascending=False)


@cached_filter
def filter_expressed_genes(
    adata: sc.AnnData,
    min_cells: int = 10,
    min_counts: float = 1.0
) -> sc.AnnData:
    """
    Filter genes by minimum expression threshold.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    min_cells : int
        Minimum number of cells with expression above threshold
    min_counts : float
        Minimum count/expression level to consider gene as expressed

    Returns
    -------
    sc.AnnData
        AnnData with low-expressed genes removed

    Notes
    -----
    For smFISH, genes need sufficient expression for detection.
    With count data, min_counts=1 means at least 1 read.
    """
    n_genes_before = adata.n_vars

    # Handle both dense and sparse matrices
    if hasattr(adata.X, 'toarray'):
        expressed_per_gene = (adata.X > min_counts).sum(axis=0).A1
    else:
        expressed_per_gene = (adata.X > min_counts).sum(axis=0)

    mask = expressed_per_gene >= min_cells
    adata = adata[:, mask].copy()

    n_genes_after = adata.n_vars
    print(f"Gene filter: {n_genes_before} → {n_genes_after} ({n_genes_before - n_genes_after} removed)")

    return adata
