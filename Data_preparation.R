# ==============================================================================
# Multiome Data Preparation for GRN inference Pipeline
# ==============================================================================
# Description: Prepares paired scRNA-seq and scATAC-seq Seurat objects into
#              a single multiome object for gene regulatory network inference
# Author: Mo Dehestani
# Date: 2025-10-10
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup and Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(
  # Input files
  rna_file = "/data/MOFA/DAN_gex.rds",
  atac_file = "/data/DAN_atac_with_motifs.rds",
  
  # Output file
  output_file = "DAN_SNCA_multiome.rds",
  
  # Cell type and condition filters
  cell_type = "DAN",        # Cell type to subset
  condition = "SNCA",       # Condition/mutation to subset
  
  # Column names in metadata
  celltype_column = "BroadCellType",
  condition_column = "Mutation",
  
  # Quality control (optional)
  min_genes = 200,          # Minimum genes per cell (set to 0 to skip)
  min_peaks = 500           # Minimum peaks per cell (set to 0 to skip)
)

# ------------------------------------------------------------------------------
# Load Required Packages
# ------------------------------------------------------------------------------

cat("Loading required packages...\n")

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

cat("✓ Packages loaded\n\n")

# ------------------------------------------------------------------------------
# Load and Subset Seurat Objects
# ------------------------------------------------------------------------------

cat("Loading Seurat objects...\n")

# Load RNA and ATAC data
rna <- readRDS(CONFIG$rna_file)
atac <- readRDS(CONFIG$atac_file)

cat(sprintf("Initial RNA: %d genes × %d cells\n", nrow(rna), ncol(rna)))
cat(sprintf("Initial ATAC: %d peaks × %d cells\n\n", nrow(atac), ncol(atac)))

# Subset by cell type
cat(sprintf("Subsetting to cell type: %s\n", CONFIG$cell_type))
Idents(rna) <- CONFIG$celltype_column
rna <- subset(rna, idents = CONFIG$cell_type)

Idents(atac) <- CONFIG$celltype_column
atac <- subset(atac, idents = CONFIG$cell_type)

# Subset by condition
cat(sprintf("Subsetting to condition: %s\n", CONFIG$condition))
Idents(rna) <- CONFIG$condition_column
rna <- subset(rna, idents = CONFIG$condition)

Idents(atac) <- CONFIG$condition_column
atac <- subset(atac, idents = CONFIG$condition)

cat(sprintf("\nAfter subsetting:\n"))
cat(sprintf("RNA: %d genes × %d cells\n", nrow(rna), ncol(rna)))
cat(sprintf("ATAC: %d peaks × %d cells\n\n", nrow(atac), ncol(atac)))

# ------------------------------------------------------------------------------
# Extract Count Matrices
# ------------------------------------------------------------------------------

cat("Extracting count matrices...\n")

countsRNA <- as.matrix(GetAssayData(rna, assay = "RNA", slot = "counts"))
countsATAC <- as.matrix(GetAssayData(atac, assay = "ATAC", slot = "counts"))

cat("✓ Count matrices extracted\n\n")

# ------------------------------------------------------------------------------
# Fix Peak Names Format
# ------------------------------------------------------------------------------

cat("Formatting peak names to chr:start-end...\n")

peaks <- rownames(countsATAC)

# Check current format
cat(sprintf("Current format example: %s\n", peaks[1]))

# Convert to chr:start-end format
if (!grepl(":", peaks[1])) {
  # Handle different separators
  if (grepl("-", peaks[1])) {
    # Format: chr1-1000-2000 → chr1:1000-2000
    peaks_split <- strsplit(peaks, "-")
    peaks_fixed <- sapply(peaks_split, function(x) {
      if(length(x) >= 3) {
        paste0(x[1], ":", x[2], "-", x[3])
      } else {
        paste(x, collapse = "-")
      }
    })
  } else if (grepl("_", peaks[1])) {
    # Format: chr1_1000_2000 → chr1:1000-2000
    peaks_split <- strsplit(peaks, "_")
    peaks_fixed <- sapply(peaks_split, function(x) {
      if(length(x) >= 3) {
        paste0(x[1], ":", x[2], "-", x[3])
      } else {
        paste(x, collapse = "_")
      }
    })
  } else {
    peaks_fixed <- peaks
    warning("Unknown peak format - keeping original names")
  }
  
  rownames(countsATAC) <- peaks_fixed
  cat(sprintf("Fixed format example: %s\n", peaks_fixed[1]))
} else {
  cat("Peak names already in correct format\n")
}

cat("✓ Peak names formatted\n\n")

# ------------------------------------------------------------------------------
# Convert Gene Symbols to Ensembl IDs
# ------------------------------------------------------------------------------

cat("Converting gene symbols to Ensembl IDs...\n")

countsRNA.df <- as.data.frame(countsRNA) %>% 
  rownames_to_column(var = "gene_symbol")

# Map symbols to Ensembl
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(countsRNA.df$gene_symbol),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "SYMBOL"
) %>%
  dplyr::filter(!is.na(ENSEMBL)) %>% 
  distinct(SYMBOL, .keep_all = TRUE)

cat(sprintf("Mapped %d genes to Ensembl IDs\n", nrow(gene_map)))

# Join and clean
countsRNA.df_mapped <- countsRNA.df %>%
  left_join(gene_map, by = c("gene_symbol" = "SYMBOL")) %>%
  dplyr::filter(!is.na(ENSEMBL)) %>%
  dplyr::select(ENSEMBL, everything(), -gene_symbol) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

cat(sprintf("Final: %d unique Ensembl genes\n", nrow(countsRNA.df_mapped)))
cat("✓ Gene mapping complete\n\n")

# ------------------------------------------------------------------------------
# Harmonize Cell Barcodes
# ------------------------------------------------------------------------------

cat("Harmonizing cell barcodes between RNA and ATAC...\n")

rna_cells <- colnames(countsRNA.df_mapped)[-1]  # Exclude ENSEMBL column
atac_cells <- colnames(countsATAC)

cat(sprintf("RNA cells: %d\n", length(rna_cells)))
cat(sprintf("ATAC cells: %d\n", length(atac_cells)))

# Define barcode transformation functions
barcode_transforms <- list(
  raw = function(x) x,
  drop_suffix_underscore = function(x) sub("_\\d+$", "", x),
  drop_suffix_twice = function(x) sub("_\\d+$", "", sub("_\\d+$", "", x)),
  drop_after_dash = function(x) sub("-\\d+_.*$", "", x),
  keep_until_dash = function(x) sub("-\\d+.*$", "-1", x),
  extract_barcode = function(x) {
    # Extract pattern: letters/numbers + dash + number
    m <- regexpr("[A-Z0-9]+-\\d+", x)
    ifelse(m > 0, regmatches(x, m), x)
  }
)

# Test all transformations
best_overlap <- 0
best_transform <- NULL
best_rna_barcodes <- NULL
best_atac_barcodes <- NULL

cat("\nTesting barcode transformations:\n")
cat(sprintf("%-30s %s\n", "Transformation", "Common Cells"))
cat(strrep("-", 50), "\n")

for (transform_name in names(barcode_transforms)) {
  transform_fn <- barcode_transforms[[transform_name]]
  
  rna_trans <- transform_fn(rna_cells)
  atac_trans <- transform_fn(atac_cells)
  common <- intersect(rna_trans, atac_trans)
  
  cat(sprintf("%-30s %d\n", transform_name, length(common)))
  
  if (length(common) > best_overlap) {
    best_overlap <- length(common)
    best_transform <- transform_name
    best_rna_barcodes <- rna_trans
    best_atac_barcodes <- atac_trans
  }
}

cat(strrep("-", 50), "\n")
cat(sprintf("\n✓ Best transformation: %s (%d common cells)\n\n", 
            best_transform, best_overlap))

# Check if we found common cells
if (best_overlap == 0) {
  cat("ERROR: No common cells found!\n\n")
  cat("First 5 RNA barcodes:\n")
  print(head(rna_cells, 5))
  cat("\nFirst 5 ATAC barcodes:\n")
  print(head(atac_cells, 5))
  cat("\nPlease check barcode formats and adjust transformations.\n")
  stop("No overlapping cells between RNA and ATAC datasets")
}

# Apply best transformation
common_cells <- intersect(best_rna_barcodes, best_atac_barcodes)

# Create mapping tables
rna_map <- data.frame(
  original = rna_cells,
  transformed = best_rna_barcodes,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(transformed %in% common_cells)

atac_map <- data.frame(
  original = atac_cells,
  transformed = best_atac_barcodes,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(transformed %in% common_cells)

# ------------------------------------------------------------------------------
# Subset to Common Cells
# ------------------------------------------------------------------------------

cat("Subsetting to common cells...\n")

# Subset RNA data
countsRNA.df_aligned <- countsRNA.df_mapped %>%
  dplyr::select(ENSEMBL, all_of(rna_map$original))

# Rename to harmonized barcodes
colnames(countsRNA.df_aligned)[-1] <- rna_map$transformed

# Subset ATAC data
countsATAC_aligned <- countsATAC[, atac_map$original]
colnames(countsATAC_aligned) <- atac_map$transformed

# Ensure same cell order
common_ordered <- sort(common_cells)

countsRNA.df_aligned <- countsRNA.df_aligned %>%
  dplyr::select(ENSEMBL, all_of(common_ordered))

countsATAC_aligned <- countsATAC_aligned[, common_ordered]

cat(sprintf("✓ Aligned: %d genes × %d cells\n", 
            nrow(countsRNA.df_aligned), length(common_ordered)))
cat(sprintf("✓ Aligned: %d peaks × %d cells\n\n", 
            nrow(countsATAC_aligned), length(common_ordered)))

# ------------------------------------------------------------------------------
# Quality Control (Optional)
# ------------------------------------------------------------------------------

if (CONFIG$min_genes > 0 || CONFIG$min_peaks > 0) {
  cat("Applying quality control filters...\n")
  
  # Calculate QC metrics
  genes_per_cell <- colSums(countsRNA.df_aligned[, -1] > 0)
  peaks_per_cell <- colSums(countsATAC_aligned > 0)
  
  # Filter cells
  cells_pass <- rep(TRUE, length(common_ordered))
  
  if (CONFIG$min_genes > 0) {
    cells_pass <- cells_pass & (genes_per_cell >= CONFIG$min_genes)
    cat(sprintf("Cells with >%d genes: %d\n", 
                CONFIG$min_genes, sum(genes_per_cell >= CONFIG$min_genes)))
  }
  
  if (CONFIG$min_peaks > 0) {
    cells_pass <- cells_pass & (peaks_per_cell >= CONFIG$min_peaks)
    cat(sprintf("Cells with >%d peaks: %d\n", 
                CONFIG$min_peaks, sum(peaks_per_cell >= CONFIG$min_peaks)))
  }
  
  if (sum(cells_pass) < length(cells_pass)) {
    cat(sprintf("\nFiltering: keeping %d / %d cells\n", 
                sum(cells_pass), length(cells_pass)))
    
    cells_keep <- common_ordered[cells_pass]
    
    countsRNA.df_aligned <- countsRNA.df_aligned %>%
      dplyr::select(ENSEMBL, all_of(cells_keep))
    
    countsATAC_aligned <- countsATAC_aligned[, cells_keep]
    common_ordered <- cells_keep
  }
  
  cat("✓ QC complete\n\n")
}

# ------------------------------------------------------------------------------
# Create Multiome Seurat Object
# ------------------------------------------------------------------------------

cat("Creating multiome Seurat object...\n")

# Extract RNA count matrix
rna_counts <- as.matrix(countsRNA.df_aligned[, -1])
rownames(rna_counts) <- countsRNA.df_aligned$ENSEMBL

# Create Seurat object with RNA
multiome <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  project = sprintf("%s_%s", CONFIG$cell_type, CONFIG$condition)
)

# Add ATAC as second assay
atac_assay <- CreateAssayObject(counts = countsATAC_aligned)
multiome[["ATAC"]] <- atac_assay

# Add metadata
multiome$cell_type <- CONFIG$cell_type
multiome$condition <- CONFIG$condition
multiome$cell_id <- colnames(multiome)

cat("✓ Multiome object created\n\n")

# ------------------------------------------------------------------------------
# Save Output
# ------------------------------------------------------------------------------

cat("Saving multiome object...\n")
saveRDS(multiome, file = CONFIG$output_file)
cat(sprintf("✓ Saved to: %s\n\n", CONFIG$output_file))

# ------------------------------------------------------------------------------
# Summary Report
# ------------------------------------------------------------------------------

cat(strrep("=", 70), "\n")
cat("MULTIOME DATA PREPARATION - SUMMARY\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("Cell type: %s\n", CONFIG$cell_type))
cat(sprintf("Condition: %s\n", CONFIG$condition))
cat(sprintf("\nFinal multiome object:\n"))
cat(sprintf("  RNA assay: %d genes\n", nrow(multiome[["RNA"]])))
cat(sprintf("  ATAC assay: %d peaks\n", nrow(multiome[["ATAC"]])))
cat(sprintf("  Cells: %d\n", ncol(multiome)))
cat(sprintf("\nBarcode harmonization: %s\n", best_transform))
cat(sprintf("Output file: %s\n", CONFIG$output_file))

# Show gene/peak statistics
cat("\nRNA expression summary:\n")
print(summary(colSums(rna_counts > 0)))

cat("\nATAC accessibility summary:\n")
print(summary(colSums(countsATAC_aligned > 0)))

cat("\n✅ Data preparation completed successfully!\n")
cat("\nNext step: Run FigR GRN inference pipeline with this multiome object.\n")