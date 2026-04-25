"""
Differential expression analysis for identifying cell type-specific markers.

This module provides functions to compute differential gene expression between
inhibitory neuron subtypes and identify genes with high specificity.

Allen VISp data:
- class/subclass/cluster columns for grouping
- Expression data is read counts (use layer='counts' or normalize first)
"""

import pandas as pd
import numpy as np
import scanpy as sc
from typing import Optional, List


def _get_expression_df(
    adata: sc.AnnData,
    use_layer: Optional[str] = None
) -> pd.DataFrame:
    """
    Extract expression matrix as DataFrame, handling sparse matrices.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    use_layer : str, optional
        AnnData layer to use (e.g., 'counts', 'normalized')
        If None, uses adata.X

    Returns
    -------
    pd.DataFrame
        Expression matrix (cells × genes)
    """
    if use_layer and use_layer in adata.layers:
        X = adata.layers[use_layer]
    else:
        X = adata.X

    # Handle sparse matrices
    if hasattr(X, 'toarray'):
        X = X.toarray()

    return pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)


def compute_cluster_means(
    adata: sc.AnnData,
    groupby: str = "subclass",
    use_layer: Optional[str] = None
) -> pd.DataFrame:
    """
    Compute mean expression for each gene in each cell type group.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    groupby : str
        Column name to group cells by (e.g., 'subclass', 'cluster')
    use_layer : str, optional
        AnnData layer to use for expression values

    Returns
    -------
    pd.DataFrame
        Mean expression matrix (genes × cell types)
    """
    expr = _get_expression_df(adata, use_layer)
    groups = adata.obs[groupby]
    mean_expr = expr.groupby(groups).mean().T

    return mean_expr


def compute_cluster_frequency(
    adata: sc.AnnData,
    groupby: str = "subclass",
    threshold: float = 0.0,
    use_layer: Optional[str] = None
) -> pd.DataFrame:
    """
    Compute detection frequency (% cells expressing) for each gene in each group.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    groupby : str
        Column name to group cells by
    threshold : float
        Expression threshold to consider gene as detected (default: 0 for counts)
    use_layer : str, optional
        AnnData layer to use for expression values

    Returns
    -------
    pd.DataFrame
        Detection frequency matrix (genes × cell types), values in [0, 100]
    """
    expr = _get_expression_df(adata, use_layer)
    groups = adata.obs[groupby]

    # Calculate detection frequency per group
    detected = (expr > threshold).astype(int)
    freq = detected.groupby(groups).mean().T * 100

    return freq


def find_differential_genes(
    adata: sc.AnnData,
    groupby: str = "subclass",
    method: str = "wilcoxon",
    n_genes: int = 100,
    use_layer: Optional[str] = None
) -> pd.DataFrame:
    """
    Identify differentially expressed genes using scanpy rank_genes_groups.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    groupby : str
        Column name defining cell groups
    method : str
        Statistical test ('wilcoxon', 't-test', 't-test_overestim_var', 'logreg')
    n_genes : int
        Number of top genes to return per group
    use_layer : str, optional
        AnnData layer to use for expression values

    Returns
    -------
    pd.DataFrame
        Differential expression results with columns:
        gene, score, pval, pval_adj, log2fc, group

    Notes
    -----
    Uses scanpy's rank_genes_groups for efficient one-vs-rest comparisons.
    For count data, consider log-normalizing first or using 'wilcoxon' method.
    """
    print(f"Running differential expression: {method}, groupby={groupby}")

    # Make a copy to avoid modifying original
    adata_de = adata.copy()

    # Run differential expression
    sc.tl.rank_genes_groups(
        adata_de,
        groupby=groupby,
        method=method,
        n_genes=n_genes,
        layer=use_layer
    )

    # Extract results
    result = adata_de.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    # Compile into DataFrame
    de_results = []
    for group in groups:
        n_genes_group = len(result['names'][group])
        group_df = pd.DataFrame({
            'gene': result['names'][group],
            'score': result['scores'][group],
            'pval': result['pvals'][group],
            'pval_adj': result['pvals_adj'][group],
            'log2fc': result['logfoldchanges'][group],
            'group': [group] * n_genes_group
        })
        de_results.append(group_df)

    de_df = pd.concat(de_results, ignore_index=True)
    print(f"Found DE genes for {len(groups)} groups")

    return de_df


def calculate_specificity_scores(
    mean_expr: pd.DataFrame,
    method: str = "gini"
) -> pd.DataFrame:
    """
    Calculate specificity scores for genes across cell types.

    Parameters
    ----------
    mean_expr : pd.DataFrame
        Mean expression matrix (genes × cell types)
    method : str
        Specificity metric: 'gini', 'tau', or 'max_ratio'

    Returns
    -------
    pd.DataFrame
        Columns: gene, specificity, top_celltype, max_expression

    Notes
    -----
    Specificity metrics:
    - gini: Gini coefficient, 0 (uniform) to 1 (cell-type specific)
    - tau: Tau index, emphasizes low expression in off-target types
    - max_ratio: Ratio of highest to second-highest expression
    """
    results = []

    for gene in mean_expr.index:
        expr = mean_expr.loc[gene].values.astype(float)

        # Skip genes with no expression
        if expr.sum() == 0:
            continue

        if method == "gini":
            sorted_expr = np.sort(expr)
            n = len(expr)
            index = np.arange(1, n + 1)
            denom = n * np.sum(sorted_expr)
            if denom > 0:
                specificity = (2 * np.sum(index * sorted_expr)) / denom - (n + 1) / n
            else:
                specificity = 0

        elif method == "tau":
            max_expr = expr.max()
            if max_expr > 0:
                specificity = np.sum(1 - expr / max_expr) / (len(expr) - 1)
            else:
                specificity = 0

        elif method == "max_ratio":
            sorted_expr = np.sort(expr)[::-1]
            if len(sorted_expr) > 1 and sorted_expr[1] > 0:
                specificity = sorted_expr[0] / sorted_expr[1]
            else:
                specificity = np.inf if sorted_expr[0] > 0 else 0

        else:
            raise ValueError(f"Unknown method: {method}. Use 'gini', 'tau', or 'max_ratio'")

        top_celltype = mean_expr.columns[expr.argmax()]

        results.append({
            'gene': gene,
            'specificity': specificity,
            'top_celltype': top_celltype,
            'max_expression': expr.max()
        })

    return pd.DataFrame(results)


def find_marker_candidates(
    adata: sc.AnnData,
    groupby: str = "subclass",
    min_mean_expression: float = 10.0,
    min_detection_freq: float = 25.0,
    min_log2fc: float = 1.0,
    max_pval_adj: float = 0.05,
    use_layer: Optional[str] = None,
    n_de_genes: int = 500
) -> pd.DataFrame:
    """
    Identify marker gene candidates based on multiple criteria.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    groupby : str
        Column defining cell type groups
    min_mean_expression : float
        Minimum mean expression (counts) in top cell type
    min_detection_freq : float
        Minimum detection frequency (%) in top cell type
    min_log2fc : float
        Minimum log2 fold-change vs other groups
    max_pval_adj : float
        Maximum adjusted p-value
    use_layer : str, optional
        AnnData layer to use
    n_de_genes : int
        Number of DE genes to compute per group

    Returns
    -------
    pd.DataFrame
        Marker candidates sorted by specificity, with columns:
        gene, celltype, mean_expression, detection_freq, log2fc, pval_adj, specificity

    Notes
    -----
    For smFISH markers, prioritizes:
    - High expression (sufficient signal)
    - High specificity (distinguishes cell types)
    - High detection frequency (consistent across cells)
    """
    print(f"Finding marker candidates for {groupby}...")
    print(f"  Thresholds: expr>={min_mean_expression}, freq>={min_detection_freq}%, "
          f"log2fc>={min_log2fc}, padj<={max_pval_adj}")

    # Compute statistics
    mean_expr = compute_cluster_means(adata, groupby=groupby, use_layer=use_layer)
    freq = compute_cluster_frequency(adata, groupby=groupby, use_layer=use_layer)
    de_results = find_differential_genes(
        adata, groupby=groupby, use_layer=use_layer, n_genes=n_de_genes
    )

    # Calculate specificity
    specificity = calculate_specificity_scores(mean_expr, method="gini")

    # Merge and filter
    markers = []
    for _, row in specificity.iterrows():
        gene = row['gene']
        celltype = row['top_celltype']

        # Check if gene is in mean_expr and freq
        if gene not in mean_expr.index or gene not in freq.index:
            continue

        # Get DE stats for this gene-celltype pair
        de_row = de_results[(de_results['gene'] == gene) & (de_results['group'] == celltype)]

        if len(de_row) == 0:
            continue

        de_row = de_row.iloc[0]

        # Get expression stats
        mean_exp = mean_expr.loc[gene, celltype]
        detection = freq.loc[gene, celltype]

        # Apply filters
        if (mean_exp >= min_mean_expression and
            detection >= min_detection_freq and
            de_row['log2fc'] >= min_log2fc and
            de_row['pval_adj'] <= max_pval_adj):

            markers.append({
                'gene': gene,
                'celltype': celltype,
                'mean_expression': mean_exp,
                'detection_freq': detection,
                'log2fc': de_row['log2fc'],
                'pval_adj': de_row['pval_adj'],
                'specificity': row['specificity']
            })

    marker_df = pd.DataFrame(markers)

    if len(marker_df) > 0:
        marker_df = marker_df.sort_values('specificity', ascending=False)

    print(f"Found {len(marker_df)} marker candidates")

    return marker_df
