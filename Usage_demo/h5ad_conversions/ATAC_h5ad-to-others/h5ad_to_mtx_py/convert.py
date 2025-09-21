#!/usr/bin/env python3
"""
H5AD to MTX Conversion Tool

This script converts ATAC-seq data from H5AD format (AnnData) to MTX format
(Matrix Market Exchange), creating sparse matrix files and metadata files for
universal compatibility across different analysis platforms.

Author: Xinze Wu
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

def h5ad_to_mtx(h5ad_file, output_dir=None):
    """
    Convert H5AD file to MTX format

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

        # Save count matrix in MTX format (cells × features - standard format)
        logger.info("Saving count matrix in MTX format")
        matrix_file = os.path.join(output_dir, "matrix.mtx")
        io.mmwrite(matrix_file, adata.X)

        # Save feature information (peaks/genes)
        logger.info("Saving feature information")
        features_df = pd.DataFrame({
            'feature_id': adata.var_names,
            'feature_name': adata.var_names,
            'feature_type': 'Peaks'
        })
        features_file = os.path.join(output_dir, "features.tsv")
        features_df.to_csv(features_file, sep='\t', header=False, index=False)

        # Save cell barcodes
        logger.info("Saving cell barcodes")
        barcodes_df = pd.DataFrame({'barcode': adata.obs_names})
        barcodes_file = os.path.join(output_dir, "barcodes.tsv")
        barcodes_df.to_csv(barcodes_file, sep='\t', header=False, index=False)

        # Save cell metadata as TSV
        if not adata.obs.empty:
            logger.info("Saving cell metadata")
            metadata_file = os.path.join(output_dir, "cell_metadata.tsv")
            adata.obs.to_csv(metadata_file, sep='\t')

        # Save feature metadata as TSV
        if not adata.var.empty:
            logger.info("Saving feature metadata")
            var_file = os.path.join(output_dir, "feature_metadata.tsv")
            adata.var.to_csv(var_file, sep='\t')

        # Save dimensionality reduction results
        save_dimensionality_reductions(adata, output_dir)

        # Create loading scripts for different platforms
        create_loading_scripts(output_dir)

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
    Save dimensionality reduction results to TSV files

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

        embedding_file = os.path.join(output_dir, f"{reduction_name}_embeddings.tsv")
        emb_df.to_csv(embedding_file, sep='\t')

def create_loading_scripts(output_dir):
    """
    Create loading scripts for different platforms

    Args:
        output_dir (str): Output directory
    """
    logger.info("Creating loading scripts for different platforms")

    # Python loading script
    python_script_content = '''# Python Script for Loading MTX Data
# Generated by h5ad_to_mtx.py

import pandas as pd
import numpy as np
from scipy import sparse, io
import scanpy as sc

def load_mtx_data(data_dir="."):
    """
    Load MTX format data back into AnnData object

    Args:
        data_dir (str): Directory containing MTX files

    Returns:
        AnnData: Loaded data object
    """
    print(f"Loading MTX data from: {data_dir}")

    # Read count matrix
    matrix = io.mmread(f"{data_dir}/matrix.mtx").tocsr()

    # Read features and barcodes
    features = pd.read_csv(f"{data_dir}/features.tsv", sep='\\t',
                          header=None, names=['feature_id', 'feature_name', 'feature_type'])
    barcodes = pd.read_csv(f"{data_dir}/barcodes.tsv", sep='\\t',
                          header=None, names=['barcode'])

    # Create AnnData object
    adata = sc.AnnData(X=matrix, dtype=np.float32)
    adata.var_names = features['feature_id'].values
    adata.obs_names = barcodes['barcode'].values

    # Load metadata if available
    try:
        cell_meta = pd.read_csv(f"{data_dir}/cell_metadata.tsv", sep='\\t', index_col=0)
        adata.obs = cell_meta.loc[adata.obs_names]
    except FileNotFoundError:
        print("No cell metadata found")

    try:
        feature_meta = pd.read_csv(f"{data_dir}/feature_metadata.tsv", sep='\\t', index_col=0)
        adata.var = feature_meta.loc[adata.var_names]
    except FileNotFoundError:
        print("No feature metadata found")

    # Load embeddings
    import glob
    embedding_files = glob.glob(f"{data_dir}/*_embeddings.tsv")
    for emb_file in embedding_files:
        reduction_name = os.path.basename(emb_file).replace('_embeddings.tsv', '')
        emb_data = pd.read_csv(emb_file, sep='\\t', index_col=0)
        common_cells = np.intersect1d(emb_data.index, adata.obs_names)
        if len(common_cells) > 0:
            adata.obsm[f'X_{reduction_name}'] = emb_data.loc[common_cells].values

    print(f"Loaded data: {adata.shape[0]} cells, {adata.shape[1]} features")
    return adata

# Example usage:
# adata = load_mtx_data()
# print(adata)
'''

    # R loading script
    r_script_content = '''# R Script for Loading MTX Data
# Generated by h5ad_to_mtx.py

library(Matrix)

load_mtx_data <- function(data_dir = ".") {
    cat("Loading MTX data from:", data_dir, "\\n")

    # Read count matrix
    matrix_path <- file.path(data_dir, "matrix.mtx")
    counts <- readMM(matrix_path)

    # Read features and barcodes
    features <- read.delim(file.path(data_dir, "features.tsv"),
                          header = FALSE, stringsAsFactors = FALSE)
    barcodes <- read.delim(file.path(data_dir, "barcodes.tsv"),
                          header = FALSE, stringsAsFactors = FALSE)

    # Set matrix row and column names
    rownames(counts) <- barcodes[, 1]
    colnames(counts) <- features[, 1]

    cat("Loaded matrix:", nrow(counts), "cells x", ncol(counts), "features\\n")

    # Load metadata if available
    metadata_list <- list()

    # Cell metadata
    cell_meta_path <- file.path(data_dir, "cell_metadata.tsv")
    if (file.exists(cell_meta_path)) {
        cell_meta <- read.delim(cell_meta_path, row.names = 1, stringsAsFactors = FALSE)
        metadata_list$cell_metadata <- cell_meta
        cat("Loaded cell metadata with", ncol(cell_meta), "columns\\n")
    }

    # Feature metadata
    feature_meta_path <- file.path(data_dir, "feature_metadata.tsv")
    if (file.exists(feature_meta_path)) {
        feature_meta <- read.delim(feature_meta_path, row.names = 1, stringsAsFactors = FALSE)
        metadata_list$feature_metadata <- feature_meta
        cat("Loaded feature metadata with", ncol(feature_meta), "columns\\n")
    }

    # Load embeddings
    embedding_files <- list.files(data_dir, pattern = "_embeddings.tsv$", full.names = TRUE)
    embeddings <- list()
    for (emb_file in embedding_files) {
        reduction_name <- gsub("_embeddings.tsv$", "", basename(emb_file))
        emb_data <- read.delim(emb_file, row.names = 1, stringsAsFactors = FALSE)
        embeddings[[reduction_name]] <- emb_data
        cat("Loaded", reduction_name, "embeddings with", ncol(emb_data), "dimensions\\n")
    }

    result <- list(
        counts = counts,
        metadata = metadata_list,
        embeddings = embeddings
    )

    cat("Data loading completed successfully!\\n")
    return(result)
}

# Example usage:
# data <- load_mtx_data()
# str(data)
'''

    # Save scripts
    python_file = os.path.join(output_dir, "load_mtx_data.py")
    with open(python_file, 'w') as f:
        f.write(python_script_content)

    r_file = os.path.join(output_dir, "load_mtx_data.R")
    with open(r_file, 'w') as f:
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

    summary_content = f"""# H5AD to MTX Conversion Summary

## Conversion Details
- **Conversion Date**: {pd.Timestamp.now()}
- **Input Format**: H5AD (AnnData)
- **Output Format**: MTX (Matrix Market Exchange)

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
- `matrix.mtx`: Sparse count matrix in Matrix Market format (cells × features)
- `features.tsv`: Feature names and types (tab-separated)
- `barcodes.tsv`: Cell barcodes/identifiers (tab-separated)
- `cell_metadata.tsv`: Cell annotations and metadata (if available)
- `feature_metadata.tsv`: Feature annotations and metadata (if available)
- `*_embeddings.tsv`: Dimensionality reduction results (tab-separated)
- `load_mtx_data.py`: Python script for loading MTX data
- `load_mtx_data.R`: R script for loading MTX data

## Platform Compatibility
The MTX format is supported by:
- **Python**: scanpy, scipy, pandas
- **R**: Matrix package, Seurat, SingleCellExperiment
- **Julia**: MatrixMarket.jl
- **MATLAB**: Matrix Market I/O functions

## Usage Examples

### Python
```python
exec(open('load_mtx_data.py').read())
adata = load_mtx_data()
```

### R
```r
source("load_mtx_data.R")
data <- load_mtx_data()
```

## Notes
- Matrix format is standard orientation (cells × features)
- Sparse matrix format is preserved for memory efficiency
- All metadata and embeddings are saved in tab-separated format
- Scripts are provided for easy data loading in Python and R
"""

    summary_file = os.path.join(output_dir, "conversion_summary.md")
    with open(summary_file, 'w') as f:
        f.write(summary_content)

def main():
    """Main function with command line interface"""
    parser = argparse.ArgumentParser(
        description='Convert H5AD (AnnData) format to MTX (Matrix Market Exchange) format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s input.h5ad                    # Convert to default output directory
  %(prog)s input.h5ad -o mtx_data        # Convert to specific directory
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
    success = h5ad_to_mtx(args.input, args.output)

    if success:
        print(f"\n✓ Conversion completed successfully!")
        print(f"Output saved to: {args.output}")
        print(f"Use the provided loading scripts to import the data")
        sys.exit(0)
    else:
        print(f"\n✗ Conversion failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()