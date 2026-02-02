# Allen Single-Cell RNA-seq Analysis for V1 L2/3 Inhibitory Neurons

This repository contains analyses of the Allen Brain Map single-cell RNA-seq dataset to identify a minimal set of RNA markers for distinguishing inhibitory cell types in Layer 2/3 of mouse primary visual cortex (V1) using serial smFISH.

## Project Goal

Identify a minimal panel of RNA species (genes) that can effectively distinguish between different inhibitory interneuron subtypes in mouse V1 Layer 2/3 for use in serial single-molecule fluorescence in situ hybridization (smFISH) experiments.

## Data Source

**Allen Brain Map - Cell Types Database**
- **URL**: https://celltypes.brain-map.org/rnaseq/mouse/v1-alm
- **Dataset**: Mouse primary visual cortex (VISp/V1) and anterior lateral motor cortex (ALM)
- **Technology**: Single-cell RNA-sequencing (Smart-seq)
- **Cell types**: Comprehensive taxonomy including excitatory and inhibitory neurons, as well as non-neuronal cells

## Scientific Background

Inhibitory interneurons in cortical Layer 2/3 comprise multiple distinct subtypes with different molecular profiles, morphologies, and functional properties. These include:
- Parvalbumin+ (Pvalb) fast-spiking interneurons
- Somatostatin+ (Sst) interneurons
- Vasoactive intestinal peptide+ (Vip) interneurons
- Additional subtypes (Lamp5, Sncg, etc.)

Serial smFISH allows for sequential detection of multiple RNA species in intact tissue, but practical constraints limit the number of probes that can be used. This analysis aims to optimize marker selection for maximal cell-type discrimination with minimal probes.

## Analysis Plan

### 1. Data Acquisition and Preprocessing
- Download Allen Brain Map scRNA-seq data for V1
- Filter for Layer 2/3 inhibitory neurons
- Quality control and normalization

### 2. Cell Type Identification
- Examine existing cell-type annotations from Allen taxonomy
- Focus on GABAergic (inhibitory) populations in L2/3
- Validate major inhibitory subtypes present in the dataset

### 3. Marker Gene Selection
- Differential expression analysis between inhibitory subtypes
- Identify genes with high specificity and expression levels
- Prioritize genes suitable for smFISH probe design:
  - High expression levels (for signal strength)
  - High specificity (minimal cross-reactivity)
  - Minimal overlap between markers

### 4. Optimization
- Determine minimal gene panel for maximal discrimination
- Consider combinatorial marker strategies
- Validate marker combinations using classifier performance

### 5. Visualization and Validation
- Generate marker expression heatmaps
- Create confusion matrices for classification accuracy
- Produce reference plots for experimental validation

## Expected Outputs

- List of candidate marker genes ranked by discriminatory power
- Recommended minimal gene panel for serial smFISH
- Expression profiles and cell-type specificity metrics
- Visualization of marker combinations and their discriminatory capacity

## Requirements

To be determined based on analysis approach (Python/R/Jupyter notebooks)

## Usage

To be added as analysis scripts are developed

## References

- Tasic, B., et al. (2018). Shared and distinct transcriptomic cell types across neocortical areas. Nature, 563(7729), 72-78.
- Yao, Z., et al. (2021). A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation. Cell, 184(12), 3222-3241.
- Allen Institute for Brain Science (2015). Allen Cell Types Database. Available from: celltypes.brain-map.org

## License

See [LICENSE](LICENSE) file for details.
