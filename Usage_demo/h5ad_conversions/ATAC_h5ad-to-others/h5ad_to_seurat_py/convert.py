#!/usr/bin/env python3
"""
H5AD to Seurat Conversion Tool

This script converts ATAC-seq data from H5AD format (AnnData) to Seurat-compatible
format, creating MTX matrix files and metadata files for seamless import into R/Seurat.

Author: ATAC Analysis Team
Version: 1.0.0
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse, io
import os
import argparse
import sys
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def validate_input_file(h5ad_file):
    """
    Validate input H5AD file exists and is readable

    Args:
        h5ad_file (str): Path to H5AD file

    Returns:
        bool: True if file is valid
    """
    if not os.path.exists(h5ad_file):
        logger.error(f"Input file does not exist: {h5ad_file}")
        return False

    if not h5ad_file.endswith('.h5ad'):
        logger.warning(f"File does not have .h5ad extension: {h5ad_file}")

    return True

def h5ad_to_seurat(h5ad_file, output_dir=None):
    """
    Convert H5AD file to Seurat-compatible format

    Args:
        h5ad_file (str): Path to input H5AD file
        output_dir (str): Output directory path (default: current directory)

    Returns:
        bool: True if conversion successful
    """
    try:
        # Use current directory if output_dir is None
        if output_dir is None:
            output_dir = os.getcwd()

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Created output directory: {output_dir}")

        # Load H5AD file
        logger.info(f"Loading H5AD file: {h5ad_file}")
        adata = sc.read_h5ad(h5ad_file)
        logger.info(f"Loaded data: {adata.shape[0]} cells, {adata.shape[1]} features")

        # Validate ATAC data format
        if not validate_atac_data(adata):
            logger.warning("Data validation issues detected, proceeding with conversion")

        # Convert sparse matrix if needed
        if not sparse.issparse(adata.X):
            logger.info("Converting to sparse matrix format")
            adata.X = sparse.csr_matrix(adata.X)

        # Save count matrix in MTX format (transposed for Seurat: features × cells)
        logger.info("Saving count matrix in MTX format")
        matrix_file = os.path.join(output_dir, "matrix.mtx")
        io.mmwrite(matrix_file, adata.X.T)

        # Save feature information (peaks/genes)
        logger.info("Saving feature information")
        features_df = pd.DataFrame({
            'gene_id': adata.var_names,
            'gene_symbol': adata.var_names,
            'feature_type': 'Peaks'
        })
        features_file = os.path.join(output_dir, "features.tsv")
        features_df.to_csv(features_file, sep='\t', header=False, index=False)

        # Save cell barcodes
        logger.info("Saving cell barcodes")
        barcodes_df = pd.DataFrame({'barcode': adata.obs_names})
        barcodes_file = os.path.join(output_dir, "barcodes.tsv")
        barcodes_df.to_csv(barcodes_file, sep='\t', header=False, index=False)

        # Save cell metadata if available
        if not adata.obs.empty:
            logger.info("Saving cell metadata")
            metadata_file = os.path.join(output_dir, "cell_metadata.csv")
            adata.obs.to_csv(metadata_file)

        # Save feature metadata if available
        if not adata.var.empty:
            logger.info("Saving feature metadata")
            var_file = os.path.join(output_dir, "feature_metadata.csv")
            adata.var.to_csv(var_file)

        # Save dimensionality reduction results
        save_dimensionality_reductions(adata, output_dir)

        # Create R import script
        create_r_import_script(output_dir)

        # Generate conversion summary
        generate_conversion_summary(adata, output_dir)

        logger.info(f"Conversion completed successfully! Files saved to: {output_dir}")
        return True

    except Exception as e:
        logger.error(f"Conversion failed: {str(e)}")
        return False

def validate_atac_data(adata):
    """
    Validate ATAC data format and quality

    Args:
        adata: AnnData object

    Returns:
        bool: True if validation passes
    """
    logger.info("Validating ATAC data format")

    # Check feature names format (should be peak coordinates)
    var_names = adata.var_names
    peak_pattern_count = sum(1 for name in var_names[:100]
                           if any(char in name for char in [':', '-', '_']))

    if peak_pattern_count / min(100, len(var_names)) < 0.5:
        logger.warning("Feature names may not be in standard peak coordinate format")

    # Check data type and range
    if sparse.issparse(adata.X):
        min_val = adata.X.min()
        max_val = adata.X.max()
    else:
        min_val = adata.X.min()
        max_val = adata.X.max()

    logger.info(f"Data range: {min_val} to {max_val}")

    if min_val < 0:
        logger.warning("Negative values detected in count matrix")

    return True

def save_dimensionality_reductions(adata, output_dir):
    """
    Save dimensionality reduction results to CSV files

    Args:
        adata: AnnData object
        output_dir (str): Output directory
    """
    for reduction_key in adata.obsm.keys():
        reduction_name = reduction_key.replace('X_', '')
        logger.info(f"Saving {reduction_name} embeddings")

        # Create DataFrame with proper column names
        n_dims = adata.obsm[reduction_key].shape[1]
        if reduction_name == 'pca':
            columns = [f'PC_{i+1}' for i in range(n_dims)]
        elif reduction_name == 'umap':
            columns = ['UMAP_1', 'UMAP_2']
        elif reduction_name == 'tsne':
            columns = ['tSNE_1', 'tSNE_2']
        else:
            columns = [f'{reduction_name}_{i+1}' for i in range(n_dims)]

        emb_df = pd.DataFrame(
            adata.obsm[reduction_key],
            index=adata.obs_names,
            columns=columns
        )

        embedding_file = os.path.join(output_dir, f"{reduction_name}_embeddings.csv")
        emb_df.to_csv(embedding_file)

def create_r_import_script(output_dir):
    """
    Create R script for importing converted data

    Args:
        output_dir (str): Output directory
    """
    logger.info("Creating R import script")

    r_script_content = '''# R Script for Importing Converted ATAC Data
# Generated by h5ad_to_seurat.py

# Load required libraries
library(Seurat)
library(Matrix)

# Function to import converted data
import_converted_data <- function(data_dir = ".") {

    cat("Importing ATAC data from:", data_dir, "\\n")

    # Read count matrix
    counts <- readMM(file.path(data_dir, "matrix.mtx"))

    # Read features and barcodes
    features <- read.delim(file.path(data_dir, "features.tsv"),
                          header = FALSE, stringsAsFactors = FALSE)
    barcodes <- read.delim(file.path(data_dir, "barcodes.tsv"),
                          header = FALSE, stringsAsFactors = FALSE)

    # Set matrix row and column names
    rownames(counts) <- features[, 1]
    colnames(counts) <- barcodes[, 1]

    # Create Seurat object
    atac_obj <- CreateSeuratObject(
        counts = counts,
        assay = "ATAC"
    )

    cat("Created Seurat object with", ncol(atac_obj), "cells and", nrow(atac_obj), "features\\n")

    # Import cell metadata if available
    metadata_file <- file.path(data_dir, "cell_metadata.csv")
    if (file.exists(metadata_file)) {
        cat("Loading cell metadata...\\n")
        cell_meta <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)
        common_cells <- intersect(rownames(cell_meta), colnames(atac_obj))
        if (length(common_cells) > 0) {
            atac_obj <- AddMetaData(atac_obj, cell_meta[common_cells, , drop = FALSE])
            cat("Added metadata for", length(common_cells), "cells\\n")
        }
    }

    # Import dimensionality reduction results
    embedding_files <- list.files(data_dir, pattern = "_embeddings.csv$", full.names = TRUE)
    for (emb_file in embedding_files) {
        reduction_name <- gsub("_embeddings.csv$", "", basename(emb_file))
        cat("Loading", reduction_name, "embeddings...\\n")

        emb_data <- read.csv(emb_file, row.names = 1, stringsAsFactors = FALSE)
        common_cells <- intersect(rownames(emb_data), colnames(atac_obj))

        if (length(common_cells) > 0) {
            # Determine key prefix
            if (reduction_name == "pca") {
                key_prefix <- "PC_"
            } else if (reduction_name == "umap") {
                key_prefix <- "UMAP_"
            } else if (reduction_name == "tsne") {
                key_prefix <- "tSNE_"
            } else {
                key_prefix <- paste0(toupper(substr(reduction_name, 1, 1)),
                                   substr(reduction_name, 2, nchar(reduction_name)), "_")
            }

            atac_obj[[reduction_name]] <- CreateDimReducObject(
                embeddings = as.matrix(emb_data[common_cells, , drop = FALSE]),
                key = key_prefix,
                assay = "ATAC"
            )
        }
    }

    cat("Data import completed successfully!\\n")
    return(atac_obj)
}

# Example usage:
# atac_data <- import_converted_data()
# print(atac_data)
# DimPlot(atac_data, reduction = "umap")  # if UMAP is available
'''

    script_file = os.path.join(output_dir, "import_seurat.R")
    with open(script_file, 'w') as f:
        f.write(r_script_content)

def generate_conversion_summary(adata, output_dir):
    """
    Generate a summary report of the conversion

    Args:
        adata: AnnData object
        output_dir (str): Output directory
    """
    logger.info("Generating conversion summary")

    # Calculate basic statistics
    if sparse.issparse(adata.X):
        sparsity = (1 - adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1])) * 100
        total_counts = adata.X.sum()
    else:
        sparsity = (np.sum(adata.X == 0) / adata.X.size) * 100
        total_counts = np.sum(adata.X)

    summary_content = f"""# H5AD to Seurat Conversion Summary

## Conversion Details
- **Conversion Date**: {pd.Timestamp.now()}
- **Input Format**: H5AD (AnnData)
- **Output Format**: Seurat-compatible (MTX + metadata)

## Data Overview
- **Cells**: {adata.shape[0]:,}
- **Features**: {adata.shape[1]:,}
- **Total Counts**: {total_counts:,}
- **Sparsity**: {sparsity:.2f}%

## Metadata Information
- **Cell Metadata Columns**: {len(adata.obs.columns)}
- **Feature Metadata Columns**: {len(adata.var.columns)}

## Dimensionality Reductions
{chr(10).join([f"- {key.replace('X_', '').upper()}: {adata.obsm[key].shape[1]} dimensions" for key in adata.obsm.keys()])}

## Output Files
- `matrix.mtx`: Sparse count matrix (features × cells)
- `features.tsv`: Feature names and types
- `barcodes.tsv`: Cell barcodes
- `cell_metadata.csv`: Cell annotations (if available)
- `feature_metadata.csv`: Feature annotations (if available)
- `*_embeddings.csv`: Dimensionality reduction results
- `import_seurat.R`: R script for data import

## Usage in R
```r
source("import_seurat.R")
atac_data <- import_converted_data()
print(atac_data)
```

## Notes
- Matrix is transposed for Seurat compatibility (features × cells)
- Sparse matrix format is preserved for memory efficiency
- All metadata and embeddings are included in the conversion
"""

    summary_file = os.path.join(output_dir, "conversion_summary.md")
    with open(summary_file, 'w') as f:
        f.write(summary_content)

def main():
    """Main function with command line interface"""
    parser = argparse.ArgumentParser(
        description='Convert H5AD (AnnData) format to Seurat-compatible format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s input.h5ad                    # Convert to default output directory
  %(prog)s input.h5ad -o seurat_data     # Convert to specific directory
  %(prog)s input.h5ad --verbose          # Enable verbose logging
        '''
    )

    parser.add_argument('input',
                       help='Input H5AD file path')
    parser.add_argument('-o', '--output',
                       default=None,
                       help='Output directory path (default: current directory)')
    parser.add_argument('--verbose',
                       action='store_true',
                       help='Enable verbose logging')

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate input
    if not validate_input_file(args.input):
        sys.exit(1)

    # Perform conversion
    success = h5ad_to_seurat(args.input, args.output)

    if success:
        print(f"\n✓ Conversion completed successfully!")
        print(f"Output saved to: {args.output}")
        print(f"Use 'source(\"import_seurat.R\")' in R to load the data")
        sys.exit(0)
    else:
        print(f"\n✗ Conversion failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()