# ==============================================================================
# Gene Regulatory Network Inference using CisTopic and FigR
# ==============================================================================
# Description: Infers transcription factor-target gene regulatory networks
#              from paired single-cell RNA-seq and ATAC-seq data
# Author: Mo Dehestani
# Date: 2025-10-10
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup and Configuration
# ------------------------------------------------------------------------------

# Set parameters
CONFIG <- list(
  # File paths
  input_file = "/data/MOFA/FigR_SNCA/DAN_SNCA_multiome.rds",
  output_dir = "figr_output",
  
  # Feature selection
  ncells = -1,              # -1 to use all cells
  nfeatures_rna = 10000,    # Number of RNA features
  nfeatures_atac = 10000,   # Number of ATAC peaks
  
  # DORC parameters
  pval_cutoff = 0.05,       # P-value cutoff for peak-gene correlations
  dorc_cutoff = 1,          # Minimum peaks per gene for DORC calling
  
  # Computational parameters
  nCores = 2,               # Number of cores for parallel processing
  kNN = 30,                 # k-nearest neighbors for smoothing
  
  # FigR parameters
  n_bg = 50,                # Background peaks for motif enrichment
  score_cut = 1             # Score cutoff for network visualization
)

# Create output directory
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Package Installation and Loading
# ------------------------------------------------------------------------------

cat("Loading required packages...\n")

# Install FigR if needed
if(!suppressMessages(require("FigR"))){
  suppressMessages(devtools::install_github("caleblareau/BuenColors"))
  suppressMessages(devtools::install_github("buenrostrolab/FigR"))
}

# Install cisTopic if needed
if(!suppressMessages(require("cisTopic")))
  suppressMessages(devtools::install_github("aertslab/cisTopic"))

# Load libraries
suppressMessages({
  library(FigR)
  library(cisTopic)
  library(Seurat)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(IRanges)  # Load early to prevent namespace conflicts
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(FNN)
  library(ComplexHeatmap)
  library(networkD3)
})

# Force correct resize function
resize <- IRanges::resize

cat("✓ Packages loaded successfully\n\n")

# ------------------------------------------------------------------------------
# Load and Prepare Data
# ------------------------------------------------------------------------------

cat("Loading multiome object...\n")
multiome <- readRDS(CONFIG$input_file)
cat(sprintf("Object class: %s\n", class(multiome)[1]))
cat(sprintf("Dimensions: %d features × %d cells\n\n", 
            nrow(multiome), ncol(multiome)))

# Normalize data
cat("Normalizing RNA and ATAC data...\n")
DefaultAssay(multiome) <- "RNA"
if(!"data" %in% Layers(multiome[["RNA"]])){
  multiome <- NormalizeData(multiome, assay = "RNA", verbose = FALSE)
}

DefaultAssay(multiome) <- "ATAC"
if(!"data" %in% Layers(multiome[["ATAC"]])){
  multiome <- RunTFIDF(multiome, assay = "ATAC", verbose = FALSE)
}

# Convert to SingleCellExperiment
RNA <- as.SingleCellExperiment(multiome, assay = "RNA")
ATAC <- as.SingleCellExperiment(multiome, assay = "ATAC")

# Ensure counts are in correct slots
if(!"counts" %in% assayNames(RNA)) counts(RNA) <- assay(RNA, "logcounts")
if(!"counts" %in% assayNames(ATAC)) counts(ATAC) <- assay(ATAC, "logcounts")

cat("✓ Data loaded and normalized\n\n")

# ------------------------------------------------------------------------------
# Feature Selection
# ------------------------------------------------------------------------------

cat("Performing feature selection...\n")

# Download TF list
tf_file <- file.path(CONFIG$output_dir, 'allTFs_hg38.txt')
if(!file.exists(tf_file)){
  download.file('https://raw.githubusercontent.com/aertslab/SCENICprotocol/master/example/allTFs_hg38.txt', 
                tf_file, quiet = TRUE)
}
tf_names <- rownames(read.table(tf_file, row.names = 1))
cat(sprintf("Loaded %d transcription factors\n", length(tf_names)))

# Subset cells if specified
if(CONFIG$ncells != -1){
  RNA <- RNA[, 1:CONFIG$ncells]
  ATAC <- ATAC[, 1:CONFIG$ncells]
}

# Select RNA features (all TFs + random non-TFs)
if(CONFIG$nfeatures_rna != -1){
  is_tf <- rownames(RNA) %in% tf_names
  index_tf <- which(is_tf)
  index_not_tf <- which(!is_tf)
  set.seed(CONFIG$nfeatures_rna)
  selected_idx <- c(index_tf, 
                    sample(index_not_tf, CONFIG$nfeatures_rna - length(index_tf)))
  RNA <- RNA[selected_idx, ]
}

# Select ATAC peaks (by accessibility variance)
if(CONFIG$nfeatures_atac != -1){
  frac_atac <- rowSums(assay(ATAC)) / ncol(ATAC)
  acc_score <- abs(0.5 - frac_atac)
  peak_mask <- rownames(ATAC) %in% names(sort(acc_score))[1:CONFIG$nfeatures_atac]
  ATAC <- ATAC[peak_mask, ]
}

cat(sprintf("Selected features: %d genes × %d cells\n", nrow(RNA), ncol(RNA)))
cat(sprintf("Selected peaks: %d peaks × %d cells\n\n", nrow(ATAC), ncol(ATAC)))

# Prepare ATAC matrix
assay(ATAC) <- as(assay(ATAC), 'sparseMatrix')
counts(ATAC) <- assay(ATAC)

# ------------------------------------------------------------------------------
# cisTopic: Topic Modeling on ATAC Data
# ------------------------------------------------------------------------------

cat("Running cisTopic for chromatin accessibility topic modeling...\n")

cistopic_file <- file.path(CONFIG$output_dir, 
                           sprintf("cisTopic_npeaks%d.rds", CONFIG$nfeatures_atac))

if(!file.exists(cistopic_file)){
  # Format peak names (chr:start-end)
  atac_df <- as.data.frame(as.matrix(counts(ATAC)))
  peak_parts <- strsplit(rownames(atac_df), "[:-]")
  
  chr <- sapply(peak_parts, `[`, 1)
  start <- sapply(peak_parts, `[`, 2)
  end <- sapply(peak_parts, `[`, 3)
  
  # Validate coordinates
  valid_idx <- !is.na(start) & !is.na(end) & 
    !is.na(as.numeric(start)) & !is.na(as.numeric(end))
  
  if(sum(!valid_idx) > 0){
    cat(sprintf("Removing %d invalid peaks\n", sum(!valid_idx)))
    atac_df <- atac_df[valid_idx, ]
    chr <- chr[valid_idx]
    start <- start[valid_idx]
    end <- end[valid_idx]
  }
  
  rownames(atac_df) <- paste0(chr, ':', start, '-', end)
  
  # Run cisTopic
  cisTopicObject <- createcisTopicObject(atac_df, project.name = 'FigR_GRN')
  cisTopicObject <- runCGSModels(cisTopicObject, 
                                 topic = 1:25,
                                 seed = 987, 
                                 nCores = CONFIG$nCores, 
                                 burnin = 90,
                                 iterations = 100, 
                                 addModels = FALSE)
  
  cisTopicObject <- selectModel(cisTopicObject, type = 'maximum')
  cisTopicObject <- runUmap(cisTopicObject, target = 'cell')
  
  topic.mat <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
  topic.mat <- t(topic.mat)
  saveRDS(topic.mat, cistopic_file)
}

cisAssign <- readRDS(cistopic_file)
cat(sprintf("✓ cisTopic complete: %d topics × %d cells\n\n", 
            nrow(cisAssign), ncol(cisAssign)))

# ------------------------------------------------------------------------------
# Cell k-NN Graph
# ------------------------------------------------------------------------------

cat("Computing cell k-nearest neighbor graph...\n")
set.seed(123)
cellkNN <- get.knn(cisAssign, k = CONFIG$kNN)$nn.index
cat(sprintf("✓ k-NN graph computed (k=%d)\n\n", CONFIG$kNN))

# ------------------------------------------------------------------------------
# Peak-Gene Correlation Analysis
# ------------------------------------------------------------------------------

cat("Computing peak-gene correlations...\n")

# Prepare RNA matrix with gene symbols
RNAmat <- as.matrix(assay(RNA))
ensembl_ids <- rownames(RNAmat)

gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Remove unmapped genes and duplicates
valid_genes <- !is.na(gene_symbols) & !duplicated(gene_symbols)
RNAmat <- RNAmat[valid_genes, ]
rownames(RNAmat) <- gene_symbols[valid_genes]

cat(sprintf("Converted to gene symbols: %d genes\n", nrow(RNAmat)))

# Prepare ATAC SummarizedExperiment
ATAC_df <- as.data.frame(as.matrix(counts(ATAC)))
peak_parts <- strsplit(rownames(ATAC_df), "[:-]")

ATAC_df$seqnames <- sapply(peak_parts, `[`, 1)
ATAC_df$start <- as.numeric(sapply(peak_parts, `[`, 2))
ATAC_df$end <- as.numeric(sapply(peak_parts, `[`, 3))

# Filter valid peaks
valid_rows <- !is.na(ATAC_df$start) & !is.na(ATAC_df$end) & 
  grepl('^chr[0-9XYM]+$', ATAC_df$seqnames) &
  ATAC_df$start < ATAC_df$end

ATAC_df <- ATAC_df[valid_rows, ]
cat(sprintf("Valid ATAC peaks: %d\n", nrow(ATAC_df)))

ATAC.se <- makeSummarizedExperimentFromDataFrame(ATAC_df)
counts(ATAC.se) <- assay(ATAC.se)
assay(ATAC.se) <- as(assay(ATAC.se), 'sparseMatrix')

# Run peak-gene correlations
ciscorr_file <- file.path(CONFIG$output_dir, 
                          sprintf("cisCorr_npeaks%d.rds", CONFIG$nfeatures_atac))

if(!file.exists(ciscorr_file)){
  cisCorr <- FigR::runGenePeakcorr(
    ATAC.se = ATAC.se,
    RNAmat = RNAmat,
    genome = "hg38",
    nCores = CONFIG$nCores,
    p.cut = NULL,
    n_bg = 250
  )
  saveRDS(cisCorr, ciscorr_file)
}

cisCorr <- readRDS(ciscorr_file)
cisCorr.filt <- cisCorr %>% dplyr::filter(pvalZ <= CONFIG$pval_cutoff)

cat(sprintf("✓ Peak-gene correlations: %d significant associations\n\n", 
            nrow(cisCorr.filt)))

# ------------------------------------------------------------------------------
# DORC Gene Identification
# ------------------------------------------------------------------------------

cat("Identifying DORC (Dynamic Regulatory Chromatin) genes...\n")

# Visualize DORC distribution
dorcGenes <- cisCorr.filt %>% 
  dorcJPlot(cutoff = CONFIG$dorc_cutoff, 
            returnGeneList = TRUE, 
            family = 'sans')

cat(sprintf("✓ Identified %d DORC genes (cutoff=%d peaks/gene)\n\n", 
            length(dorcGenes), CONFIG$dorc_cutoff))

if(length(dorcGenes) < 30){
  warning("Low number of DORC genes. Consider lowering dorc_cutoff or pval_cutoff.")
}

# ------------------------------------------------------------------------------
# DORC Score Calculation and Smoothing
# ------------------------------------------------------------------------------

cat("Computing DORC scores and smoothing across k-NN...\n")

dorcMat <- getDORCScores(ATAC.se, 
                         dorcTab = cisCorr.filt, 
                         geneList = dorcGenes, 
                         nCores = CONFIG$nCores)

rownames(cellkNN) <- colnames(dorcMat)
dorcMat.s <- smoothScoresNN(NNmat = cellkNN[, 1:20], 
                            mat = dorcMat, 
                            nCores = CONFIG$nCores)

# Smooth RNA expression
rownames(cellkNN) <- colnames(RNAmat)
RNAmat.s <- smoothScoresNN(NNmat = cellkNN[, 1:20], 
                           mat = RNAmat, 
                           nCores = CONFIG$nCores)

cat(sprintf("✓ Smoothed matrices: %d DORC genes × %d cells\n\n", 
            nrow(dorcMat.s), ncol(dorcMat.s)))

# ------------------------------------------------------------------------------
# FigR: Gene Regulatory Network Inference
# ------------------------------------------------------------------------------

cat("Running FigR to infer transcription factor-target relationships...\n")

dorcK <- max(3, min(10, nrow(dorcMat.s) - 2))

figr_file <- file.path(CONFIG$output_dir, "figR_GRN_results.rds")

if(!file.exists(figr_file)){
  figR.d <- runFigRGRN(
    ATAC.se = ATAC.se,
    dorcTab = cisCorr.filt,
    genome = "hg38",
    dorcMat = dorcMat.s,
    rnaMat = RNAmat.s,
    dorcK = dorcK,
    n_bg = CONFIG$n_bg,
    nCores = CONFIG$nCores
  )
  
  saveRDS(figR.d, figr_file)
  write.csv(figR.d, 
            file.path(CONFIG$output_dir, "figR_GRN_results.csv"), 
            row.names = FALSE)
}

figR.d <- readRDS(figr_file)

cat(sprintf("✓ FigR complete: %d TF-target interactions identified\n\n", 
            nrow(figR.d)))

# ------------------------------------------------------------------------------
# Visualization and Analysis
# ------------------------------------------------------------------------------

cat("Generating visualizations...\n")

# 1. TF-DORC Enrichment Scatter Plot
library(ggrastr)
library(BuenColors)

p1 <- figR.d %>%
  ggplot(aes(Corr.log10P, Enrichment.log10P, color = Score)) +
  ggrastr::geom_point_rast(size = 0.5, shape = 16) +
  theme_classic() +
  scale_color_gradientn(colours = jdb_palette("solar_extra"), 
                        limits = c(-3, 3), 
                        oob = scales::squish) +
  labs(title = "TF-DORC Regulatory Relationships",
       x = "Peak-Gene Correlation (-log10 P)",
       y = "Motif Enrichment (-log10 P)")

ggsave(file.path(CONFIG$output_dir, "figR_scatter.png"), 
       p1, width = 8, height = 6, dpi = 300)

# 2. Rank Driver TFs
drivers <- rankDrivers(figR.d, 
                       rankBy = "meanScore", 
                       interactive = FALSE)

write.csv(drivers, 
          file.path(CONFIG$output_dir, "ranked_driver_TFs.csv"), 
          row.names = FALSE)

cat(sprintf("Top 10 Driver TFs:\n"))
print(head(drivers, 10))

# 3. Regulatory Network Heatmap
png(file.path(CONFIG$output_dir, "figR_heatmap.png"), 
    width = 1200, height = 1000, res = 150)

heatmap <- plotfigRHeatmap(
  figR.d = figR.d,
  score.cut = CONFIG$score_cut,
  TFs = unique(figR.d$Motif),
  show_row_dend = FALSE
)

draw(heatmap, newpage = TRUE)
dev.off()

# 4. Interactive Network
top_tfs <- names(sort(table(figR.d$Motif), decreasing = TRUE)[1:10])

d3_network <- plotfigRNetwork(
  figR.d,
  score.cut = CONFIG$score_cut,
  weight.edges = TRUE,
  TFs = top_tfs
)

htmlwidgets::saveWidget(d3_network, 
                        file.path(CONFIG$output_dir, "figR_network.html"),
                        selfcontained = TRUE)

cat("\n✓ All visualizations saved\n\n")

# ------------------------------------------------------------------------------
# Summary Report
# ------------------------------------------------------------------------------

cat("=" %R% 78, "\n")
cat("FIGR GENE REGULATORY NETWORK ANALYSIS - SUMMARY\n")
cat("=" %R% 78, "\n\n")

cat(sprintf("Input cells: %d\n", ncol(multiome)))
cat(sprintf("RNA features analyzed: %d\n", nrow(RNA)))
cat(sprintf("ATAC peaks analyzed: %d\n", nrow(ATAC)))
cat(sprintf("Significant peak-gene correlations: %d\n", nrow(cisCorr.filt)))
cat(sprintf("DORC genes identified: %d\n", length(dorcGenes)))
cat(sprintf("TF-target interactions: %d\n", nrow(figR.d)))
cat(sprintf("Unique TFs: %d\n", length(unique(figR.d$Motif))))
cat(sprintf("Unique target genes: %d\n", length(unique(figR.d$Gene))))
cat(sprintf("\nOutput directory: %s\n", CONFIG$output_dir))

cat("\n✅ Pipeline completed successfully!\n")
