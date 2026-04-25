"""
Visualization utilities for scRNA-seq marker analysis.

This module provides plotting functions for exploring cell types,
marker expression, and panel performance.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from typing import List, Optional, Tuple


def plot_celltype_distribution(
    adata: sc.AnnData,
    groupby: str = "subclass",
    figsize: Tuple[int, int] = (10, 6)
) -> plt.Figure:
    """
    Plot distribution of cells across cell types.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    groupby : str
        Column to group by
    figsize : Tuple[int, int]
        Figure size

    Returns
    -------
    plt.Figure
        Matplotlib figure
    """
    counts = adata.obs[groupby].value_counts()

    fig, ax = plt.subplots(figsize=figsize)
    counts.plot(kind='barh', ax=ax)
    ax.set_xlabel('Number of cells')
    ax.set_ylabel('Cell type')
    ax.set_title(f'Cell type distribution ({len(adata)} cells)')

    plt.tight_layout()
    return fig


def plot_marker_heatmap(
    adata: sc.AnnData,
    genes: List[str],
    groupby: str = "subclass",
    layer: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 8),
    cmap: str = "RdYlBu_r"
) -> plt.Figure:
    """
    Plot heatmap of marker expression across cell types.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    genes : List[str]
        Genes to plot
    groupby : str
        Column defining cell type groups
    layer : str, optional
        Layer to use for expression values
    figsize : Tuple[int, int]
        Figure size
    cmap : str
        Colormap name

    Returns
    -------
    plt.Figure
        Matplotlib figure
    """
    # Get expression data
    if layer:
        expr = pd.DataFrame(
            adata[:, genes].layers[layer],
            index=adata.obs_names,
            columns=genes
        )
    else:
        expr = pd.DataFrame(
            adata[:, genes].X,
            index=adata.obs_names,
            columns=genes
        )

    # Compute mean expression per cell type
    groups = adata.obs[groupby]
    mean_expr = expr.groupby(groups).mean().T

    # Plot heatmap
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        mean_expr,
        cmap=cmap,
        cbar_kws={'label': 'Mean expression'},
        ax=ax
    )
    ax.set_title('Marker gene expression by cell type')
    ax.set_xlabel('Cell type')
    ax.set_ylabel('Gene')

    plt.tight_layout()
    return fig


def plot_marker_dotplot(
    adata: sc.AnnData,
    genes: List[str],
    groupby: str = "subclass",
    layer: Optional[str] = None,
    figsize: Tuple[int, int] = (12, 6)
) -> plt.Figure:
    """
    Create dotplot showing marker expression and detection frequency.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    genes : List[str]
        Genes to plot
    groupby : str
        Column defining cell type groups
    layer : str, optional
        Layer to use for expression values
    figsize : Tuple[int, int]
        Figure size

    Returns
    -------
    plt.Figure
        Matplotlib figure with dotplot
    """
    # Use scanpy's dotplot
    fig = sc.pl.dotplot(
        adata,
        genes,
        groupby=groupby,
        layer=layer,
        return_fig=True,
        figsize=figsize
    )

    return fig


def plot_confusion_matrix(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    labels: Optional[List[str]] = None,
    figsize: Tuple[int, int] = (10, 8)
) -> plt.Figure:
    """
    Plot confusion matrix for classification results.

    Parameters
    ----------
    y_true : np.ndarray
        True labels
    y_pred : np.ndarray
        Predicted labels
    labels : List[str], optional
        Label names
    figsize : Tuple[int, int]
        Figure size

    Returns
    -------
    plt.Figure
        Matplotlib figure
    """
    from sklearn.metrics import confusion_matrix

    cm = confusion_matrix(y_true, y_pred, labels=labels)

    # Normalize to percentages
    cm_norm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis] * 100

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        cm_norm,
        annot=True,
        fmt='.1f',
        cmap='Blues',
        xticklabels=labels,
        yticklabels=labels,
        ax=ax,
        cbar_kws={'label': 'Percentage (%)'}
    )
    ax.set_xlabel('Predicted')
    ax.set_ylabel('True')
    ax.set_title('Classification confusion matrix')

    plt.tight_layout()
    return fig


def plot_separability_heatmap(
    separability_df: pd.DataFrame,
    figsize: Tuple[int, int] = (10, 8)
) -> plt.Figure:
    """
    Plot pairwise cell type separability heatmap.

    Parameters
    ----------
    separability_df : pd.DataFrame
        Pairwise separability matrix
    figsize : Tuple[int, int]
        Figure size

    Returns
    -------
    plt.Figure
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        separability_df,
        annot=True,
        fmt='.2f',
        cmap='RdYlGn',
        vmin=0.5,
        vmax=1.0,
        center=0.75,
        ax=ax,
        cbar_kws={'label': 'Classification accuracy'}
    )
    ax.set_title('Pairwise cell type separability')

    plt.tight_layout()
    return fig


def plot_gene_expression_violin(
    adata: sc.AnnData,
    genes: List[str],
    groupby: str = "subclass",
    layer: Optional[str] = None,
    figsize: Tuple[int, int] = (12, 4)
) -> plt.Figure:
    """
    Plot violin plots for gene expression across cell types.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    genes : List[str]
        Genes to plot
    groupby : str
        Column defining cell type groups
    layer : str, optional
        Layer to use for expression values
    figsize : Tuple[int, int]
        Figure size

    Returns
    -------
    plt.Figure
        Matplotlib figure
    """
    fig = sc.pl.violin(
        adata,
        genes,
        groupby=groupby,
        layer=layer,
        return_fig=True
    )

    return fig


def plot_umap_with_markers(
    adata: sc.AnnData,
    genes: List[str],
    layer: Optional[str] = None,
    ncols: int = 3,
    figsize: Optional[Tuple[int, int]] = None
) -> plt.Figure:
    """
    Plot UMAP colored by marker gene expression.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object (must have UMAP coordinates)
    genes : List[str]
        Genes to plot
    layer : str, optional
        Layer to use for expression values
    ncols : int
        Number of columns in subplot grid
    figsize : Tuple[int, int], optional
        Figure size

    Returns
    -------
    plt.Figure
        Matplotlib figure
    """
    if 'X_umap' not in adata.obsm:
        print("Computing UMAP coordinates...")
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

    fig = sc.pl.umap(
        adata,
        color=genes,
        layer=layer,
        ncols=ncols,
        return_fig=True,
        frameon=False
    )

    return fig

"""
Python port of the Allen Institute's `varibow` colormap generator.

Original R implementation lives in scrattch.io/R/misc.R:
    https://github.com/AllenInstitute/scrattch.io/blob/master/R/misc.R

Generates a palette by cycling hue (like rainbow) while ALSO cycling
saturation through {0.55, 0.7, 0.85, 1.0} and value through {1.0, 0.8, 0.6}.
The coprime cycle lengths (4 and 3) make adjacent colors easy to tell apart
even with hundreds of cell types.
"""

import colorsys


def varibow(n_colors):
    """Return a list of n_colors hex strings (e.g. '#FF0000')."""
    sats = [0.55, 0.70, 0.85, 1.00]
    vals = [1.00, 0.80, 0.60]

    # R's rainbow() walks hue from 0 to (n-1)/n, so hue_i = i / n
    hexes = []
    for i in range(n_colors):
        h = i / n_colors
        s = sats[i % len(sats)]
        v = vals[i % len(vals)]
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        hexes.append("#{:02X}{:02X}{:02X}".format(
            round(r * 255), round(g * 255), round(b * 255)
        ))
    return hexes


# if __name__ == "__main__":
#     # Quick demo: print 10 colors and render a swatch PNG
#     palette = varibow(60)
#     for i, c in enumerate(palette):
#         print(f"{i:3d}  {c}")

#     try:
#         import matplotlib.pyplot as plt
#         from matplotlib.patches import Rectangle

#         n = 60
#         palette = varibow(n)
#         fig, ax = plt.subplots(figsize=(10, 1.2))
#         for i, c in enumerate(palette):
#             ax.add_patch(Rectangle((i, 0), 1, 1, facecolor=c, edgecolor="none"))
#         ax.set_xlim(0, n); ax.set_ylim(0, 1)
#         ax.set_xticks([]); ax.set_yticks([])
#         ax.set_title(f"varibow({n})")
#         plt.tight_layout()
#         # plt.savefig("/home/claude/varibow_swatch.png", dpi=120)
#         print("\nSaved swatch to varibow_swatch.png")
#     except ImportError:
#         pass