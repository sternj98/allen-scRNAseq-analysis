# AnnData Object Reference — `adata_l23_inh`

**Source**: Allen Brain Map VISp SMART-seq (2018-06-14)  
**Loaded via**: `data_loader.load_allen_filtered()`  
**Filter**: GABAergic cells in superficial layers (`L1`, `L1-L2/3`, `L1-L4`, `L1-L6`, `L2/3`, `L2/3-L4`, `L4`)

```
AnnData: 3215 cells × 45768 genes
```

---

## Core arrays

| Attribute | Shape | dtype | Contents |
|---|---|---|---|
| `.X` | 3215 × 45768 | float32 | Raw exon counts (copy of `layers['counts']`) |
| `.layers['counts']` | 3215 × 45768 | float32 | Raw exon counts — use this for analysis to preserve the unmodified values if `.X` is normalized |

---

## `.obs` — per-cell metadata (3215 rows)

Indexed by Allen `sample_name` (e.g. `F1S4_160108_001_A01`).

### Cell-type annotations
| Column | Description | Example values |
|---|---|---|
| `class` | Broad cell class | `GABAergic` |
| `subclass` | Inhibitory subclass | `Vip`, `Lamp5`, `Pvalb`, `Sst`, `Sncg`, `Serpinf1`, `Meis2` |
| `cluster` | Fine-grained cluster | `Vip Chat Htr1f`, `Pvalb Tpbg`, … (50 total, 36 with ≥10 cells) |
| `core_intermediate_call` | Assignment confidence | `Core`, `Intermediate` |
| `confusion_score` | Classifier ambiguity (lower = better) | float |
| `cluster_correlation` | Correlation to cluster centroid | float |

### Spatial / anatomical
| Column | Description | Example values |
|---|---|---|
| `brain_region` | Gross region | `VISp` |
| `brain_subregion` | Cortical layer | `L2/3`, `L4`, `L1-L6`, … |
| `brain_hemisphere` | Left / right | `right` |

### Donor / sample metadata
`sample_id`, `sample_type`, `organism`, `donor`, `sex`, `age_days`, `eye_condition`, `genotype`, `driver_lines`, `reporter_lines`

### Library / sequencing QC
`facs_date`, `facs_container`, `facs_sort_criteria`, `rna_amplification_set`, `library_prep_set`, `library_prep_avg_size_bp`, `seq_name`, `seq_tube`, `seq_batch`

### Expression QC metrics
| Column | Description |
|---|---|
| `total_reads` | Total mapped reads |
| `percent_exon_reads` | % reads mapping to exons |
| `percent_intron_reads` | % reads mapping to introns |
| `percent_intergenic_reads` | % intergenic reads |
| `percent_rrna_reads` | % ribosomal RNA reads |
| `percent_mt_exon_reads` | % mitochondrial exon reads |
| `percent_reads_unique` | % uniquely mapped reads |
| `percent_synth_reads` | % synthetic spike-in reads |
| `percent_ecoli_reads` | % E. coli contamination reads |
| `percent_aligned_reads_total` | Total alignment rate |
| `complexity_cg` | Library complexity (CG content) |
| `genes_detected_cpm_criterion` | Genes detected at CPM threshold |
| `genes_detected_fpkm_criterion` | Genes detected at FPKM threshold |
| `tdt_cpm` | TdTomato reporter CPM |
| `gfp_cpm` | GFP reporter CPM |

---

## `.var` — per-gene metadata (45768 rows)

Indexed by `gene_entrez_id` (as string, e.g. `"71661"`). This matches the expression matrix row index in the Allen CSV files.

| Column | Description | Example |
|---|---|---|
| `gene_symbol` | HGNC/MGI gene symbol | `Gad1`, `Pvalb`, `Vip` |
| `gene_id` | Allen internal gene ID | `500717483` |
| `chromosome` | Chromosome | `2`, `X` |
| `gene_entrez_id` | NCBI Entrez ID (= index) | `71661` |
| `gene_name` | Full gene name | `glutamate decarboxylase 1` |

**Looking up a gene by symbol:**
```python
symbol_to_id = pd.Series(adata_l23_inh.var_names, index=adata_l23_inh.var['gene_symbol'])
entrez_id = symbol_to_id['Pvalb']
adata_l23_inh[:, entrez_id]  # subset to Pvalb column
```

**Gene ID → symbol mapping** (used throughout analysis notebooks):
```python
id_to_symbol = adata_l23_inh.var['gene_symbol']   # Series: entrez_id → symbol
```

---

## `.obsm` — per-cell embeddings

| Key | Shape | Contents |
|---|---|---|
| `X_umap` | 3215 × 2 | UMAP coordinates computed on all GABAergic cells (6125), sliced to this subset. Loaded from `{DATA_DIR}/cache/umap_gaba.csv`. |

---

## `.uns` — unstructured metadata

| Key | Contents |
|---|---|
| `_cache_metadata` | Dict with `function`, `timestamp`, `params`, `n_input_cells`, `n_output_cells` — written by `cell_filter.cached_filter` when a result is loaded from cache. |

---

## Subclass composition

| Subclass | N cells |
|---|---|
| Vip | 1327 |
| Lamp5 | 981 |
| Pvalb | 513 |
| Sst | 315 |
| Sncg | 66 |
| Serpinf1 | 7 |
| Meis2 | 6 |

36 clusters have ≥10 cells and are used for mapping analysis.

---

## Common indexing patterns

```python
# Subset to a subclass
adata_vip = adata_l23_inh[adata_l23_inh.obs['subclass'] == 'Vip']

# Get raw counts for panel genes
panel_col_idx = [adata_l23_inh.var_names.get_loc(gid) for gid in panel_ids]
X_panel = adata_l23_inh.layers['counts'][:, panel_col_idx]

# UMAP coordinates aligned to a cell subset
cell_positions = adata_l23_inh.obs_names.get_indexer(obs_subset.index)
umap_subset = adata_l23_inh.obsm['X_umap'][cell_positions]
```
