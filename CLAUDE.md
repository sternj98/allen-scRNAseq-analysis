# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a computational biology project analyzing Allen Brain Map single-cell RNA-seq data to identify minimal marker gene panels for distinguishing inhibitory neuron subtypes in mouse V1 Layer 2/3 using serial smFISH.

**Key Scientific Goal**: Identify the smallest set of RNA species that can reliably distinguish between different GABAergic interneuron subtypes in cortical Layer 2/3 for serial single-molecule FISH experiments.

**Data Source**: Allen Brain Map Cell Types Database (https://celltypes.brain-map.org/rnaseq/mouse/v1-alm)

**Target Cell Types**: Inhibitory interneurons in Layer 2/3, including Pvalb+, Sst+, Vip+, Lamp5+, and other GABAergic subtypes.

## Project Structure

This is a new repository. The expected structure will include:
- Analysis notebooks or scripts for scRNA-seq data processing
- Data directories (likely gitignored) for Allen dataset downloads
- Output directories for figures, tables, and marker gene lists
- Documentation of marker selection criteria and validation results

## Analysis Workflow

### Expected Analysis Steps
1. **Data Download**: Obtain Allen Brain Map scRNA-seq metadata and expression matrices for V1
2. **Cell Filtering**: Subset to Layer 2/3 inhibitory neurons based on Allen taxonomy annotations
3. **Differential Expression**: Identify genes that distinguish between inhibitory subtypes
4. **Marker Selection**: Prioritize genes by expression level, specificity, and suitability for smFISH
5. **Optimization**: Determine minimal gene panel that maximizes cell-type discrimination
6. **Validation**: Test classifier performance with selected marker combinations

### Typical Python/R Environments
This project will likely use:
- **Python**: scanpy, pandas, numpy, matplotlib/seaborn for scRNA-seq analysis
- **Jupyter notebooks** for exploratory analysis and visualization

## Domain-Specific Considerations

### Single-Cell RNA-seq Analysis
- Allen data uses Smart-seq2 technology (full-length transcript sequencing)
- Cell-type annotations are pre-computed by Allen Institute
- Expression values may be in TPM, CPM, or log-normalized formats
- Focus on genes with sufficient expression for smFISH detection (typically >1-2 CPM/TPM)

### smFISH Marker Requirements
When selecting marker genes, prioritize:
- **High expression**: Genes with robust signal in target cell types
- **Specificity**: Minimal expression in off-target populations
- **Non-overlapping**: Markers that distinguish complementary cell groups
- **Probe-suitable**: Genes with sufficient mRNA length and GC content for probe design
- **Validated**: Genes with known expression patterns in literature when possible

### Allen Data Format
- Allen provides pre-processed expression matrices and comprehensive metadata
- Cell-type annotations follow hierarchical taxonomy (class > subclass > cluster)
- Layer assignments and regional annotations are included in metadata
- May need to map between gene symbols and Ensembl IDs

## Common Commands

*To be added as analysis pipeline develops*

## Architecture Notes

*To be populated as code structure emerges*

When this repository develops analysis code, update this section with:
- Key modules/scripts and their purposes
- Data flow between analysis steps
- Helper functions for marker selection and validation
- Visualization pipelines
