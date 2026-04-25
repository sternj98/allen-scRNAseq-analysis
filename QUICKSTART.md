# Quick Start Guide

## Setup

1. **Install dependencies**
```bash
pip install -r requirements.txt
```

2. **Download Allen data**
   - Visit https://celltypes.brain-map.org/rnaseq/mouse/v1-alm
   - Download metadata and expression matrix
   - Save to `data/raw/` directory

## Option 1: Interactive Analysis (Recommended)

Open the Jupyter notebook for step-by-step analysis with visualizations:

```bash
jupyter notebook notebooks/01_marker_discovery.ipynb
```

Follow the notebook cells to:
- Load and explore data
- Filter to Layer 2/3 inhibitory neurons
- Identify marker candidates
- Optimize marker panel
- Validate results

## Option 2: Command-Line Script

Run the complete pipeline automatically:

```bash
python run_analysis.py \
    --metadata data/raw/allen_metadata.csv \
    --expression data/raw/allen_expression.csv \
    --output outputs/markers/ \
    --max-markers 10 \
    --min-accuracy 0.95
```

**Options:**
- `--metadata`: Path to Allen metadata CSV
- `--expression`: Path to expression matrix CSV
- `--output`: Output directory for results
- `--class-col`: Column name for cell class (default: 'class')
- `--layer-col`: Column name for layer (default: 'layer')
- `--subclass-col`: Column name for subclass (default: 'subclass')
- `--max-markers`: Maximum markers to select (default: 10)
- `--min-accuracy`: Target accuracy (default: 0.95)

## Option 3: Python API

Use the modules programmatically:

```python
import sys
sys.path.append('src')

from data_loader import load_allen_data, create_anndata
from cell_filter import filter_layer_23_inhibitory
from diff_expression import find_marker_candidates
from marker_selection import greedy_marker_selection

# Load and filter data
metadata, expression = load_allen_data('data/raw/allen_metadata.csv',
                                       'data/raw/allen_expression.csv')
adata = create_anndata(metadata, expression)
adata_l23 = filter_layer_23_inhibitory(adata)

# Find markers
markers = find_marker_candidates(adata_l23, groupby='subclass')
optimal = greedy_marker_selection(adata_l23,
                                  markers['gene'].tolist(),
                                  max_markers=10)
print(f"Optimal markers: {optimal}")
```

## Expected Outputs

After running the analysis, you'll find:

**Marker Lists:**
- `outputs/markers/optimal_marker_panel.csv` - Final selected markers
- `outputs/markers/all_marker_candidates.csv` - All candidates
- `outputs/markers/pairwise_separability.csv` - Cell type discrimination

**Figures** (from notebook):
- Cell type distributions
- Marker expression heatmaps
- Separability matrices
- UMAP visualizations

## Troubleshooting

**"No cells matched filters"**
- Check that column names match your data
- Use `--class-col`, `--layer-col`, `--subclass-col` to specify correct columns
- Verify layer annotations include '2/3' or 'L2/3'

**Low marker panel accuracy**
- Increase `--max-markers` to allow more genes
- Decrease `--min-accuracy` target
- Check that sufficient cells exist per subclass

**Memory issues**
- Filter genes more aggressively (increase `min_expression`)
- Process subsets of cell types separately
- Use sparse matrix format in scanpy

## Next Steps

1. Review marker expression patterns in the notebook
2. Examine pairwise separability to identify problematic cell type pairs
3. Consider tissue-specific expression patterns from literature
4. Design smFISH probes for selected markers
5. Validate markers experimentally

## Resources

- [Allen Cell Types Database](https://celltypes.brain-map.org/)
- [Scanpy Documentation](https://scanpy.readthedocs.io/)
- [smFISH Probe Design Guidelines](https://www.biosearchtech.com/support/education/stellaris-rna-fish)
