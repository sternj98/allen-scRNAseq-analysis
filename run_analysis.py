#!/usr/bin/env python
"""
Command-line script to run the complete marker discovery analysis.

Usage:
    python run_analysis.py --metadata data/raw/allen_metadata.csv \
                          --expression data/raw/allen_expression.csv \
                          --output outputs/markers/
"""

import argparse
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

import pandas as pd
import scanpy as sc
import data_loader
import cell_filter
import diff_expression
import marker_selection
import visualization


def main():
    parser = argparse.ArgumentParser(
        description='Run Allen scRNA-seq marker discovery analysis'
    )
    parser.add_argument(
        '--metadata',
        type=str,
        required=True,
        help='Path to Allen metadata CSV file'
    )
    parser.add_argument(
        '--expression',
        type=str,
        required=True,
        help='Path to Allen expression matrix CSV file'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='outputs/markers/',
        help='Output directory for results'
    )
    parser.add_argument(
        '--class-col',
        type=str,
        default='class',
        help='Column name for cell class'
    )
    parser.add_argument(
        '--layer-col',
        type=str,
        default='layer',
        help='Column name for layer annotation'
    )
    parser.add_argument(
        '--subclass-col',
        type=str,
        default='subclass',
        help='Column name for cell subclass'
    )
    parser.add_argument(
        '--max-markers',
        type=int,
        default=10,
        help='Maximum number of markers to select'
    )
    parser.add_argument(
        '--min-accuracy',
        type=float,
        default=0.95,
        help='Target classification accuracy'
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("Allen Brain Map Marker Discovery Pipeline")
    print("=" * 80)

    # 1. Load data
    print("\n[1/6] Loading Allen data...")
    metadata, expression = data_loader.load_allen_data(
        args.metadata,
        args.expression
    )
    adata = data_loader.create_anndata(metadata, expression)
    print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # 2. Filter to Layer 2/3 inhibitory neurons
    print("\n[2/6] Filtering to Layer 2/3 inhibitory neurons...")
    adata_l23 = cell_filter.filter_layer_23_inhibitory(
        adata,
        class_col=args.class_col,
        layer_col=args.layer_col
    )

    # Filter genes
    adata_l23 = cell_filter.filter_expressed_genes(
        adata_l23,
        min_cells=10,
        min_expression=1.0
    )

    # Show cell type distribution
    summary = cell_filter.get_cell_type_summary(adata_l23, groupby=args.subclass_col)
    print("\nCell type distribution:")
    print(summary)

    # 3. Find marker candidates
    print("\n[3/6] Identifying marker candidates...")
    marker_candidates = diff_expression.find_marker_candidates(
        adata_l23,
        groupby=args.subclass_col,
        min_mean_expression=2.0,
        min_detection_freq=30.0,
        min_log2fc=1.0,
        max_pval_adj=0.05
    )

    # Save all candidates
    candidates_path = output_dir / 'all_marker_candidates.csv'
    marker_candidates.to_csv(candidates_path, index=False)
    print(f"Saved {len(marker_candidates)} candidates to {candidates_path}")

    # 4. Get top markers per cell type
    print("\n[4/6] Selecting markers per cell type...")
    celltype_markers = marker_selection.optimize_panel_per_celltype(
        adata_l23,
        marker_candidates,
        groupby=args.subclass_col,
        n_markers_per_type=2
    )

    # Compile candidate pool
    candidate_pool = []
    for markers in celltype_markers.values():
        candidate_pool.extend(markers)
    candidate_pool = list(set(candidate_pool))
    print(f"\nCandidate pool: {len(candidate_pool)} genes")

    # 5. Optimize marker panel
    print("\n[5/6] Optimizing marker panel...")
    optimal_markers = marker_selection.greedy_marker_selection(
        adata_l23,
        candidate_genes=candidate_pool,
        groupby=args.subclass_col,
        max_markers=args.max_markers,
        min_accuracy=args.min_accuracy
    )

    print(f"\nOptimal marker panel ({len(optimal_markers)} genes):")
    print(", ".join(optimal_markers))

    # 6. Validate and export
    print("\n[6/6] Validating and exporting results...")

    # Score panel
    final_accuracy = marker_selection.score_marker_panel(
        adata_l23,
        genes=optimal_markers,
        groupby=args.subclass_col,
        cv=10
    )
    print(f"Final panel accuracy (10-fold CV): {final_accuracy:.4f}")

    # Compute pairwise separability
    separability = marker_selection.calculate_pairwise_separability(
        adata_l23,
        genes=optimal_markers,
        groupby=args.subclass_col
    )

    # Export results
    panel_path = output_dir / 'optimal_marker_panel.csv'
    marker_selection.export_marker_panel(
        genes=optimal_markers,
        output_path=str(panel_path),
        adata=adata_l23,
        marker_df=marker_candidates,
        groupby=args.subclass_col
    )

    sep_path = output_dir / 'pairwise_separability.csv'
    separability.to_csv(sep_path)
    print(f"Saved pairwise separability to {sep_path}")

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == '__main__':
    main()
