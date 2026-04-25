# Allen Brain Map scRNA-seq Marker Discovery

Computational pipeline for identifying minimal marker gene panels to distinguish Layer 2/3 GABAergic interneuron subtypes for serial smFISH experiments.

## Project Structure

```
├── data/
│   ├── raw/              # Downloaded Allen data (gitignored)
│   └── processed/        # Processed datasets (gitignored)
├── notebooks/
│   └── 01_marker_discovery.ipynb  # Main analysis workflow
├── src/
│   ├── data_loader.py    # Data download and loading utilities
│   ├── cell_filter.py    # Cell filtering functions
│   ├── diff_expression.py  # Differential expression analysis
│   ├── marker_selection.py  # Marker panel optimization
│   └── visualization.py  # Plotting utilities
├── outputs/
│   ├── figures/          # Generated plots
│   └── markers/          # Marker gene lists
└── requirements.txt      # Python dependencies
```

## Installation

1. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Data Download

Download Allen Brain Map scRNA-seq data from:
https://celltypes.brain-map.org/rnaseq/mouse/v1-alm

Required files:
- **Cell metadata**: Contains cell type annotations, layer assignments, and cluster IDs
- **Expression matrix**: Gene expression values (TPM or counts)

Save downloaded files to:
- `data/raw/allen_metadata.csv`
- `data/raw/allen_expression.csv`

## Usage

### Quick Start

Run the main analysis notebook:
```bash
jupyter notebook notebooks/01_marker_discovery.ipynb
```

### Module Usage

```python
import sys
sys.path.append('src')

import data_loader
import cell_filter
import diff_expression
import marker_selection

# Load data
metadata, expression = data_loader.load_allen_data(
    'data/raw/allen_metadata.csv',
    'data/raw/allen_expression.csv'
)

# Create AnnData object
adata = data_loader.create_anndata(metadata, expression)

# Filter to Layer 2/3 inhibitory neurons
adata_l23 = cell_filter.filter_layer_23_inhibitory(
    adata,
    class_col='class',
    layer_col='layer'
)

# Find marker candidates
markers = diff_expression.find_marker_candidates(
    adata_l23,
    groupby='subclass',
    min_mean_expression=2.0,
    min_detection_freq=30.0
)

# Optimize marker panel
optimal_markers = marker_selection.greedy_marker_selection(
    adata_l23,
    candidate_genes=markers['gene'].tolist(),
    max_markers=10,
    min_accuracy=0.95
)
```

## Analysis Workflow

### 1. Data Loading and QC
- Load Allen Brain Map scRNA-seq data
- Create AnnData object for analysis
- Apply basic quality control filters

### 2. Cell Filtering
- Subset to Layer 2/3 neurons
- Filter to inhibitory (GABAergic) cells
- Remove low-quality cells and rare clusters

### 3. Differential Expression
- Compute cell type-specific expression statistics
- Identify differentially expressed genes
- Calculate specificity scores

### 4. Marker Selection
- Filter candidates by smFISH criteria:
  - High expression (>2 TPM)
  - High detection frequency (>30% cells)
  - Statistical significance (FDR < 0.05)
- Optimize minimal marker panel using greedy selection
- Validate classification performance

### 5. Validation
- Cross-validated accuracy assessment
- Pairwise cell type separability analysis
- Expression pattern visualization

## Output Files

### Marker Gene Lists
- `outputs/markers/optimal_marker_panel.csv`: Final selected markers with statistics
- `outputs/markers/all_marker_candidates.csv`: All candidate markers passing filters
- `outputs/markers/pairwise_separability.csv`: Cell type discrimination matrix

### Figures
- `outputs/figures/final_marker_heatmap.png`: Expression heatmap of selected markers
- Additional plots generated during exploratory analysis

## Key Functions

### Data Loading
- `load_allen_data()`: Load metadata and expression from CSV
- `create_anndata()`: Create AnnData object
- `filter_by_quality()`: Apply QC filters

### Cell Filtering
- `filter_layer_23_inhibitory()`: Subset to Layer 2/3 GABAergic neurons
- `filter_by_subclass()`: Select specific inhibitory subtypes
- `filter_expressed_genes()`: Remove low-expressed genes

### Differential Expression
- `find_differential_genes()`: One-vs-rest differential expression
- `calculate_specificity_scores()`: Compute gene specificity metrics
- `find_marker_candidates()`: Identify smFISH-suitable markers

### Marker Selection
- `greedy_marker_selection()`: Forward feature selection
- `exhaustive_panel_search()`: Test all combinations (small panels)
- `score_marker_panel()`: Evaluate classification accuracy
- `calculate_pairwise_separability()`: Compute cell type discrimination

## Notes

### smFISH Marker Criteria
Selected markers are optimized for smFISH detection:
- **Expression level**: ≥2 TPM for reliable signal
- **Detection frequency**: ≥30% cells in target type
- **Specificity**: High fold-change vs other types
- **Statistical support**: FDR < 0.05

### Allen Data Format
- **Smart-seq2 protocol**: Full-length transcript sequencing
- **Hierarchical taxonomy**: class → subclass → cluster
- **Pre-computed annotations**: Cell types and layers already assigned
- **Expression units**: Typically TPM (transcripts per million)

## References

- Allen Institute Cell Types Database: https://celltypes.brain-map.org/
- Scanpy: https://scanpy.readthedocs.io/
- Tasic et al. (2018) Nature: Shared and distinct transcriptomic cell types across neocortical areas

## Contact

For questions about this analysis pipeline, see CLAUDE.md for project details.
