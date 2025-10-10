
# Single-cell Multiomics GRN Inference Pipeline

Infers transcription factor (TF)-target gene regulatory networks from paired single-cell RNA-seq and ATAC-seq data using the FigR algorithm.

## Overview

This pipeline processes multiome (paired scRNA-seq + scATAC-seq) data to identify gene regulatory relationships by:

1. **Topic modeling** chromatin accessibility patterns (cisTopic)
2. **Correlating** ATAC peaks with gene expression 
3. **Identifying** DORC (Dynamic Open Regulatory Chromatin) genes
4. **Inferring** TF-target relationships via motif enrichment in correlated peaks

## Requirements

### R Dependencies
```r
# Bioconductor packages
BSgenome.Hsapiens.UCSC.hg38
org.Hs.eg.db
SingleCellExperiment
ComplexHeatmap

# GitHub packages
FigR (buenrostrolab/FigR)
cisTopic (aertslab/cisTopic)
BuenColors (caleblareau/BuenColors)

# CRAN packages
Seurat, dplyr, ggplot2, Matrix, IRanges
FNN, networkD3, ggrastr
```

### Input Data
- **Seurat multiome object** (`.rds`) containing:
  - `RNA` assay: scRNA-seq counts
  - `ATAC` assay: scATAC-seq peak counts
  - Both assays normalized (or script will normalize)

## Quick Start

```r
# 1. Edit configuration
CONFIG <- list(
  input_file = "path/to/multiome.rds",
  output_dir = "figr_output",
  nfeatures_rna = 10000,
  nfeatures_atac = 10000,
  nCores = 4
)

# 2. Run pipeline
source("figr_grn_pipeline.R")
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nfeatures_rna` | 10000 | Number of genes to analyze (keeps all TFs) |
| `nfeatures_atac` | 10000 | Number of ATAC peaks to analyze |
| `pval_cutoff` | 0.05 | P-value threshold for peak-gene correlations |
| `dorc_cutoff` | 1 | Minimum peaks per gene for DORC calling |
| `kNN` | 30 | k-neighbors for expression smoothing |
| `score_cut` | 1 | Minimum FigR score for visualization |

**Tip**: If you get too few DORC genes, try:
- Lowering `dorc_cutoff` to 0.5
- Increasing `pval_cutoff` to 0.1
- Setting `nfeatures_atac = -1` (use all peaks)

## Output Files

```
figr_output/
├── cisTopic_npeaks10000.rds          # Topic modeling results
├── cisCorr_npeaks10000.rds           # Peak-gene correlations
├── figR_GRN_results.rds              # Main GRN results (R object)
├── figR_GRN_results.csv              # GRN results (table)
├── ranked_driver_TFs.csv             # Top regulatory TFs
├── figR_scatter.png                  # TF-DORC enrichment plot
├── figR_heatmap.png                  # TF-target heatmap
└── figR_network.html                 # Interactive network visualization
```

## Key Results

### GRN Table (`figR_GRN_results.csv`)
Columns include:
- `Motif`: Transcription factor
- `Gene`: Target gene
- `Peak`: Associated ATAC peak
- `Score`: FigR regulatory score
- `Cor`: Peak-gene correlation
- `Enrichment.log10P`: Motif enrichment significance

### Interpretation
- **High Score**: Strong evidence for TF→gene regulation
- **Positive Cor**: Peak accessibility correlates with gene expression
- **High Enrichment**: TF motif significantly enriched in peak

## Method Details

**FigR Score** = Peak-gene correlation × Motif enrichment

The pipeline:
1. Models chromatin accessibility topics across cells (cisTopic)
2. Computes peak-gene correlations within cis-regulatory domains (±250kb)
3. Identifies DORC genes (genes with multiple correlated peaks)
4. Tests TF motif enrichment in DORC peaks
5. Scores TF-target relationships by correlation × enrichment

## Troubleshooting

### "Error in resize()..."
- **Cause**: Namespace conflict between IRanges and cisTopic
- **Fix**: Script loads `IRanges` early and sets `resize <- IRanges::resize`

### Too few DORC genes (<30)
- Lower `dorc_cutoff` from 1 to 0.5
- Increase `pval_cutoff` from 0.05 to 0.1
- Use more peaks: `nfeatures_atac = -1`

### Long runtime
- Reduce `nfeatures_atac` and `nfeatures_rna`
- Increase `nCores` for parallel processing
- Use subset of cells: `ncells = 1000`

## Citation

If you use this pipeline, please cite the tools on which this pipeline relys:

**FigR**: Kartha et al. (2022) "Functional inference of gene regulation using single-cell multi-omics" *Cell Genomics*

**cisTopic**: Bravo González-Blas et al. (2019) "cisTopic: cis-regulatory topic modeling on single-cell ATAC-seq data" *Nature Methods*

## License

MIT License - See LICENSE file for details

## Author

Mo Dehestani - Oct 2025

## Repository

https://github.com/yourusername/figr-grn-pipeline
