"""
Marker gene selection and panel optimization.

This module provides functions to optimize minimal marker gene panels for
distinguishing cell types, with focus on smFISH experimental constraints.

Key functions:
- score_marker_panel: Evaluate classification accuracy of a gene panel
- greedy_marker_selection: Forward feature selection for optimal panel
- calculate_pairwise_separability: Identify difficult cell type pairs
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from itertools import combinations
import scanpy as sc


def _get_expression_matrix(
    adata: sc.AnnData,
    genes: List[str],
    use_layer: Optional[str] = None
) -> np.ndarray:
    """
    Extract expression matrix for specified genes.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData
    genes : List[str]
        Genes to extract
    use_layer : str, optional
        AnnData layer to use

    Returns
    -------
    np.ndarray
        Expression matrix (cells × genes)
    """
    # Filter to genes that exist in adata
    valid_genes = [g for g in genes if g in adata.var_names]
    if len(valid_genes) < len(genes):
        missing = set(genes) - set(valid_genes)
        print(f"Warning: {len(missing)} genes not found: {list(missing)[:5]}...")

    if len(valid_genes) == 0:
        raise ValueError("No valid genes found in adata")

    if use_layer and use_layer in adata.layers:
        X = adata[:, valid_genes].layers[use_layer]
    else:
        X = adata[:, valid_genes].X

    if hasattr(X, 'toarray'):
        X = X.toarray()

    return X


def score_marker_panel(
    adata: sc.AnnData,
    genes: List[str],
    groupby: str = "subclass",
    method: str = "random_forest",
    cv: int = 5,
    use_layer: Optional[str] = None
) -> float:
    """
    Score a marker gene panel by classification accuracy.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    genes : List[str]
        Marker genes to evaluate
    groupby : str
        Column defining cell type labels
    method : str
        Classifier: 'random_forest' or 'logistic'
    cv : int
        Number of cross-validation folds
    use_layer : str, optional
        AnnData layer to use for expression

    Returns
    -------
    float
        Mean cross-validated classification accuracy
    """
    X = _get_expression_matrix(adata, genes, use_layer)
    y = adata.obs[groupby].values

    # Use stratified CV to handle class imbalance
    cv_splitter = StratifiedKFold(n_splits=cv, shuffle=True, random_state=42)

    if method == "random_forest":
        clf = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            n_jobs=-1
        )
    else:
        from sklearn.linear_model import LogisticRegression
        clf = LogisticRegression(
            max_iter=1000,
            random_state=42,
            n_jobs=-1
        )

    scores = cross_val_score(clf, X, y, cv=cv_splitter, scoring='accuracy')
    return scores.mean()


def greedy_marker_selection(
    adata: sc.AnnData,
    candidate_genes: List[str],
    groupby: str = "subclass",
    max_markers: int = 10,
    min_accuracy: float = 0.95,
    use_layer: Optional[str] = None,
    verbose: bool = True
) -> List[str]:
    """
    Greedily select minimal marker panel using forward feature selection.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    candidate_genes : List[str]
        Pool of candidate marker genes
    groupby : str
        Column defining cell type labels
    max_markers : int
        Maximum number of markers to select
    min_accuracy : float
        Target classification accuracy (stopping criterion)
    use_layer : str, optional
        AnnData layer to use
    verbose : bool
        Print progress

    Returns
    -------
    List[str]
        Selected marker genes in order of addition
    """
    if verbose:
        print(f"Greedy selection: max={max_markers}, target={min_accuracy}")
        print(f"Candidate pool: {len(candidate_genes)} genes")

    # Filter to valid genes
    valid_candidates = [g for g in candidate_genes if g in adata.var_names]
    if len(valid_candidates) < len(candidate_genes):
        print(f"  {len(candidate_genes) - len(valid_candidates)} candidates not in adata")

    selected = []
    remaining = list(valid_candidates)
    best_score = 0.0

    for i in range(max_markers):
        best_gene = None
        best_new_score = best_score

        # Try adding each remaining gene
        for gene in remaining:
            test_panel = selected + [gene]
            try:
                score = score_marker_panel(
                    adata, test_panel, groupby=groupby, use_layer=use_layer, cv=5
                )
                if score > best_new_score:
                    best_new_score = score
                    best_gene = gene
            except Exception as e:
                continue

        if best_gene is None:
            if verbose:
                print(f"  Step {i+1}: No improvement found")
            break

        selected.append(best_gene)
        remaining.remove(best_gene)
        best_score = best_new_score

        if verbose:
            print(f"  Step {i+1}: +{best_gene} → accuracy={best_score:.4f}")

        if best_score >= min_accuracy:
            if verbose:
                print(f"  Reached target accuracy {min_accuracy}")
            break

    return selected


def exhaustive_panel_search(
    adata: sc.AnnData,
    candidate_genes: List[str],
    groupby: str = "subclass",
    panel_size: int = 5,
    use_layer: Optional[str] = None,
    top_n: int = 10
) -> pd.DataFrame:
    """
    Exhaustively search all possible marker panels of given size.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    candidate_genes : List[str]
        Pool of candidate genes (keep small, <20)
    groupby : str
        Column defining cell type labels
    panel_size : int
        Number of genes per panel
    use_layer : str, optional
        AnnData layer to use
    top_n : int
        Number of top panels to return

    Returns
    -------
    pd.DataFrame
        Top panels with columns: genes, accuracy
    """
    valid_candidates = [g for g in candidate_genes if g in adata.var_names]
    n_combos = len(list(combinations(valid_candidates, panel_size)))

    print(f"Testing {n_combos} panels of size {panel_size}")

    if n_combos > 5000:
        print(f"WARNING: {n_combos} combinations - this may take a while!")

    results = []
    for i, genes in enumerate(combinations(valid_candidates, panel_size)):
        if (i + 1) % 100 == 0:
            print(f"  Progress: {i+1}/{n_combos}")

        try:
            score = score_marker_panel(
                adata, list(genes), groupby=groupby, use_layer=use_layer
            )
            results.append({'genes': list(genes), 'accuracy': score})
        except Exception:
            continue

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('accuracy', ascending=False).head(top_n)

    return results_df.reset_index(drop=True)


def calculate_pairwise_separability(
    adata: sc.AnnData,
    genes: List[str],
    groupby: str = "subclass",
    use_layer: Optional[str] = None,
    min_cells: int = 10
) -> pd.DataFrame:
    """
    Calculate pairwise separability between all cell type pairs.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    genes : List[str]
        Marker genes to evaluate
    groupby : str
        Column defining cell type labels
    use_layer : str, optional
        AnnData layer to use
    min_cells : int
        Minimum cells per type for comparison

    Returns
    -------
    pd.DataFrame
        Symmetric matrix of pairwise classification accuracy
    """
    X = _get_expression_matrix(adata, genes, use_layer)

    # Get cell types with sufficient cells
    type_counts = adata.obs[groupby].value_counts()
    valid_types = type_counts[type_counts >= min_cells].index.tolist()

    print(f"Computing pairwise separability for {len(valid_types)} cell types")

    separability = pd.DataFrame(
        np.ones((len(valid_types), len(valid_types))),
        index=valid_types,
        columns=valid_types
    )

    for i, ct1 in enumerate(valid_types):
        for j, ct2 in enumerate(valid_types):
            if i >= j:
                continue

            # Subset to these two cell types
            mask = adata.obs[groupby].isin([ct1, ct2])
            X_pair = X[mask.values]
            y_pair = adata.obs.loc[mask, groupby].values

            # Skip if too few cells
            if len(y_pair) < 20:
                separability.loc[ct1, ct2] = np.nan
                separability.loc[ct2, ct1] = np.nan
                continue

            try:
                clf = RandomForestClassifier(n_estimators=50, random_state=42, n_jobs=-1)
                score = cross_val_score(clf, X_pair, y_pair, cv=3, scoring='accuracy').mean()
                separability.loc[ct1, ct2] = score
                separability.loc[ct2, ct1] = score
            except Exception:
                separability.loc[ct1, ct2] = np.nan
                separability.loc[ct2, ct1] = np.nan

    return separability


def optimize_panel_per_celltype(
    adata: sc.AnnData,
    marker_df: pd.DataFrame,
    groupby: str = "subclass",
    n_markers_per_type: int = 2
) -> Dict[str, List[str]]:
    """
    Select top markers for each cell type independently.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object
    marker_df : pd.DataFrame
        Marker candidates with 'gene' and 'celltype' columns
    groupby : str
        Column defining cell type labels
    n_markers_per_type : int
        Number of markers per cell type

    Returns
    -------
    Dict[str, List[str]]
        Mapping of cell type to list of top markers
    """
    panel = {}

    for celltype in adata.obs[groupby].unique():
        ct_markers = marker_df[marker_df['celltype'] == celltype].copy()

        if len(ct_markers) == 0:
            panel[celltype] = []
            print(f"{celltype}: (no markers found)")
            continue

        ct_markers = ct_markers.sort_values(
            ['specificity', 'mean_expression'],
            ascending=[False, False]
        )

        top_markers = ct_markers.head(n_markers_per_type)['gene'].tolist()
        panel[celltype] = top_markers

        print(f"{celltype}: {', '.join(top_markers) if top_markers else '(none)'}")

    return panel


def get_combined_marker_pool(
    panel_dict: Dict[str, List[str]]
) -> List[str]:
    """
    Combine markers from all cell types into a single unique list.

    Parameters
    ----------
    panel_dict : Dict[str, List[str]]
        Output from optimize_panel_per_celltype

    Returns
    -------
    List[str]
        Unique list of all markers
    """
    all_markers = []
    for markers in panel_dict.values():
        all_markers.extend(markers)
    return list(set(all_markers))


def export_marker_panel(
    genes: List[str],
    output_path: str,
    adata: sc.AnnData,
    marker_df: pd.DataFrame,
    groupby: str = "subclass"
) -> pd.DataFrame:
    """
    Export marker panel with statistics to CSV.

    Parameters
    ----------
    genes : List[str]
        Selected marker genes
    output_path : str
        Path to save CSV
    adata : sc.AnnData
        Input AnnData
    marker_df : pd.DataFrame
        Marker candidate statistics
    groupby : str
        Column defining cell types

    Returns
    -------
    pd.DataFrame
        Exported data
    """
    export_data = []

    for gene in genes:
        gene_info = marker_df[marker_df['gene'] == gene]

        if len(gene_info) > 0:
            info = gene_info.iloc[0]
            export_data.append({
                'gene': gene,
                'top_celltype': info['celltype'],
                'mean_expression': info['mean_expression'],
                'detection_freq': info['detection_freq'],
                'log2fc': info['log2fc'],
                'pval_adj': info['pval_adj'],
                'specificity': info['specificity']
            })
        else:
            export_data.append({'gene': gene})

    export_df = pd.DataFrame(export_data)

    # Score the panel
    try:
        accuracy = score_marker_panel(adata, genes, groupby=groupby)
        print(f"Panel accuracy: {accuracy:.4f}")
    except Exception as e:
        print(f"Could not score panel: {e}")

    export_df.to_csv(output_path, index=False)
    print(f"Exported to {output_path}")

    return export_df
