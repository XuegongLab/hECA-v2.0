#!/usr/bin/env python3
"""
H5AD to CSV Conversion Tool (Fixed Version v2.0)

This script converts ATAC-seq data from H5AD format (AnnData) to CSV format,
creating files compatible with various analysis platforms and tools that prefer
tabular data formats.

Fixed: f-string syntax error and added gzip compression support

Author: ATAC Analysis Team
Version: 2.0.0
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import os
import argparse
import sys
import logging
import gzip
from pathlib import Path

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

def estimate_memory_usage(adata):
    """
    Estimate memory usage for CSV conversion

    Args:
        adata: AnnData object

    Returns:
        dict: Memory estimates
    """
    n_cells, n_features = adata.shape

    # Estimate dense matrix size (worst case)
    dense_size_gb = (n_cells * n_features * 8) / (1024**3)  # 8 bytes per float64

    # Estimate CSV file size (approximate)
    avg_digits_per_value = 6  # Average digits per number
    csv_size_gb = (n_cells * n_features * avg_digits_per_value) / (1024**3)

    # Estimate compressed size (roughly 70-90% compression for sparse data)
    if sparse.issparse(adata.X):
        sparsity = (1 - adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1]))
        compression_ratio = 0.1 + (sparsity * 0.8)  # Better compression for sparser data
    else:
        compression_ratio = 0.3  # Standard compression for dense data

    compressed_size_gb = csv_size_gb * compression_ratio

    return {
        'dense_matrix_gb': dense_size_gb,
        'csv_file_gb': csv_size_gb,
        'compressed_csv_gb': compressed_size_gb,
        'sparsity': sparsity if sparse.issparse(adata.X) else 0.0
    }

def h5ad_to_csv(h5ad_file, output_dir=None, compress=True, max_size_mb=1000):
    """
    Convert H5AD file to CSV format with optional compression

    Args:
        h5ad_file (str): Path to input H5AD file
        output_dir (str): Output directory path
        compress (bool): Whether to compress CSV files with gzip
        max_size_mb (int): Maximum uncompressed file size in MB before splitting

    Returns:
        bool: True if conversion successful
    """
    try:
        # Set default output directory in current working directory
        if output_dir is None:
            base_name = os.path.splitext(os.path.basename(h5ad_file))[0]
            output_dir = os.path.join(os.getcwd(), base_name + "_csv")

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Created output directory: {output_dir}")

        # Load H5AD file
        logger.info(f"Loading H5AD file: {h5ad_file}")
        adata = sc.read_h5ad(h5ad_file)
        logger.info(f"Loaded data: {adata.shape[0]} cells, {adata.shape[1]} features")

        # Estimate memory usage
        memory_est = estimate_memory_usage(adata)
        logger.info(f"Memory estimates:")
        logger.info(f"  Dense matrix: {memory_est['dense_matrix_gb']:.2f} GB")
        logger.info(f"  CSV file: {memory_est['csv_file_gb']:.2f} GB")
        if compress:
            logger.info(f"  Compressed CSV: {memory_est['compressed_csv_gb']:.2f} GB")
        logger.info(f"  Data sparsity: {memory_est['sparsity']:.2%}")

        # Check if file will be too large
        if memory_est['csv_file_gb'] * 1024 > max_size_mb:
            logger.warning(f"Estimated CSV size ({memory_est['csv_file_gb']:.2f} GB) exceeds limit ({max_size_mb} MB)")
            logger.info("Consider using compression or splitting options")

        # Convert sparse matrix to dense for CSV export
        logger.info("Converting sparse matrix to dense format")
        if sparse.issparse(adata.X):
            matrix = adata.X.toarray()
        else:
            matrix = adata.X.copy()

        # Create main count matrix DataFrame
        logger.info("Creating count matrix DataFrame")
        count_df = pd.DataFrame(
            matrix,
            index=adata.obs_names,
            columns=adata.var_names
        )

        # Save count matrix
        count_file = "count_matrix.csv"
        if compress:
            count_file += ".gz"

        count_path = os.path.join(output_dir, count_file)
        logger.info(f"Saving count matrix to: {count_path}")

        if compress:
            count_df.to_csv(count_path, compression='gzip')
        else:
            count_df.to_csv(count_path)

        # Save cell metadata
        if not adata.obs.empty:
            logger.info("Saving cell metadata")
            cell_meta_file = "cell_metadata.csv"
            if compress:
                cell_meta_file += ".gz"

            cell_meta_path = os.path.join(output_dir, cell_meta_file)

            if compress:
                adata.obs.to_csv(cell_meta_path, compression='gzip')
            else:
                adata.obs.to_csv(cell_meta_path)
        else:
            logger.info("No cell metadata to save")

        # Save feature metadata
        if not adata.var.empty:
            logger.info("Saving feature metadata")
            feature_meta_file = "feature_metadata.csv"
            if compress:
                feature_meta_file += ".gz"

            feature_meta_path = os.path.join(output_dir, feature_meta_file)

            if compress:
                adata.var.to_csv(feature_meta_path, compression='gzip')
            else:
                adata.var.to_csv(feature_meta_path)
        else:
            logger.info("No feature metadata to save")

        # Save embeddings
        if adata.obsm.keys():
            logger.info("Saving dimensionality reduction embeddings")
            for key in adata.obsm.keys():
                embedding_name = key.replace('X_', '')
                embedding_file = f"embedding_{embedding_name}.csv"
                if compress:
                    embedding_file += ".gz"

                embedding_path = os.path.join(output_dir, embedding_file)

                embedding_df = pd.DataFrame(
                    adata.obsm[key],
                    index=adata.obs_names,
                    columns=[f"{embedding_name}_{i+1}" for i in range(adata.obsm[key].shape[1])]
                )

                if compress:
                    embedding_df.to_csv(embedding_path, compression='gzip')
                else:
                    embedding_df.to_csv(embedding_path)

                logger.info(f"Saved {embedding_name} embedding ({adata.obsm[key].shape[1]} dimensions)")

        # Create loading scripts
        create_loading_scripts(output_dir, compress)

        # Generate conversion summary
        generate_conversion_summary(adata, output_dir, compress)

        logger.info("Conversion completed successfully!")
        return True

    except Exception as e:
        logger.error(f"Conversion failed: {str(e)}")
        logger.exception("Full error traceback:")
        return False

def create_loading_scripts(output_dir, compress):
    """
    Create scripts for loading CSV data back into analysis environments

    Args:
        output_dir (str): Output directory path
        compress (bool): Whether files are compressed
    """
    logger.info("Creating loading scripts")

    # File extensions
    ext = ".csv.gz" if compress else ".csv"

    # Python loading script
    python_script_content = f'''# Python Loading Script for CSV Data (v2.0)
# Generated by h5ad_to_csv.py v2.0

import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path

def load_csv_data(data_dir="."):
    """
    Load CSV data back into Python for analysis

    Args:
        data_dir (str): Directory containing CSV files

    Returns:
        dict: Dictionary containing loaded data
    """
    data_dir = Path(data_dir)
    data = {{}}

    # Load count matrix
    count_file = data_dir / "count_matrix{ext}"
    if count_file.exists():
        print("Loading count matrix...")
        if "{ext}" == ".csv.gz":
            data['count_matrix'] = pd.read_csv(count_file, compression='gzip', index_col=0)
        else:
            data['count_matrix'] = pd.read_csv(count_file, index_col=0)
        print(f"Count matrix shape: {{data['count_matrix'].shape}}")

    # Load cell metadata
    cell_meta_file = data_dir / "cell_metadata{ext}"
    if cell_meta_file.exists():
        print("Loading cell metadata...")
        if "{ext}" == ".csv.gz":
            data['cell_metadata'] = pd.read_csv(cell_meta_file, compression='gzip', index_col=0)
        else:
            data['cell_metadata'] = pd.read_csv(cell_meta_file, index_col=0)
        print(f"Cell metadata shape: {{data['cell_metadata'].shape}}")

    # Load feature metadata
    feature_meta_file = data_dir / "feature_metadata{ext}"
    if feature_meta_file.exists():
        print("Loading feature metadata...")
        if "{ext}" == ".csv.gz":
            data['feature_metadata'] = pd.read_csv(feature_meta_file, compression='gzip', index_col=0)
        else:
            data['feature_metadata'] = pd.read_csv(feature_meta_file, index_col=0)
        print(f"Feature metadata shape: {{data['feature_metadata'].shape}}")

    # Load embeddings
    embedding_files = list(data_dir.glob(f"embedding_*{ext}"))
    if embedding_files:
        print("Loading embeddings...")
        data['embeddings'] = {{}}
        for emb_file in embedding_files:
            emb_name = emb_file.stem.replace("embedding_", "")
            if emb_name.endswith(".csv"):  # Remove .csv if present
                emb_name = emb_name[:-4]

            if "{ext}" == ".csv.gz":
                emb_data = pd.read_csv(emb_file, compression='gzip', index_col=0)
            else:
                emb_data = pd.read_csv(emb_file, index_col=0)

            data['embeddings'][emb_name] = emb_data
            print(f"Loaded {{emb_name}} embedding: {{emb_data.shape}}")

    return data

def create_anndata_object(data):
    """
    Create AnnData object from loaded CSV data

    Args:
        data (dict): Dictionary from load_csv_data()

    Returns:
        AnnData: Reconstructed AnnData object
    """
    if 'count_matrix' not in data:
        raise ValueError("Count matrix is required")

    # Create AnnData object
    adata = sc.AnnData(X=data['count_matrix'].values)
    adata.obs_names = data['count_matrix'].index
    adata.var_names = data['count_matrix'].columns

    # Add cell metadata
    if 'cell_metadata' in data:
        # Align cell metadata with count matrix
        common_cells = adata.obs_names.intersection(data['cell_metadata'].index)
        if len(common_cells) > 0:
            adata.obs = data['cell_metadata'].loc[common_cells]

    # Add feature metadata
    if 'feature_metadata' in data:
        # Align feature metadata with count matrix
        common_features = adata.var_names.intersection(data['feature_metadata'].index)
        if len(common_features) > 0:
            adata.var = data['feature_metadata'].loc[common_features]

    # Add embeddings
    if 'embeddings' in data:
        for emb_name, emb_data in data['embeddings'].items():
            # Align embeddings with count matrix
            common_cells = adata.obs_names.intersection(emb_data.index)
            if len(common_cells) > 0:
                adata.obsm[f'X_{{emb_name}}'] = emb_data.loc[common_cells].values

    print(f"Created AnnData object: {{adata.shape[0]}} cells, {{adata.shape[1]}} features")
    return adata

# Example usage:
# data = load_csv_data(".")
# adata = create_anndata_object(data)
'''

    # R loading script - fix the f-string backslash issue
    newline = "\\n"  # Define newline outside f-string
    r_script_content = f'''# R Loading Script for CSV Data (v2.0)
# Generated by h5ad_to_csv.py v2.0

library(data.table)  # Fast CSV reading
library(Seurat)

# Function to load CSV data
load_csv_data <- function(data_dir = ".") {{
    cat("Loading CSV data from:", data_dir, "{newline}")

    data <- list()

    # Load count matrix
    count_file <- file.path(data_dir, "count_matrix{ext}")
    if (file.exists(count_file)) {{
        cat("Loading count matrix...{newline}")
        data$count_matrix <- as.data.frame(fread(count_file))
        rownames(data$count_matrix) <- data$count_matrix[[1]]
        data$count_matrix[[1]] <- NULL
        cat("Count matrix shape:", dim(data$count_matrix), "{newline}")
    }}

    # Load cell metadata
    cell_meta_file <- file.path(data_dir, "cell_metadata{ext}")
    if (file.exists(cell_meta_file)) {{
        cat("Loading cell metadata...{newline}")
        data$cell_metadata <- as.data.frame(fread(cell_meta_file))
        rownames(data$cell_metadata) <- data$cell_metadata[[1]]
        data$cell_metadata[[1]] <- NULL
        cat("Cell metadata shape:", dim(data$cell_metadata), "{newline}")
    }}

    # Load feature metadata
    feature_meta_file <- file.path(data_dir, "feature_metadata{ext}")
    if (file.exists(feature_meta_file)) {{
        cat("Loading feature metadata...{newline}")
        data$feature_metadata <- as.data.frame(fread(feature_meta_file))
        rownames(data$feature_metadata) <- data$feature_metadata[[1]]
        data$feature_metadata[[1]] <- NULL
        cat("Feature metadata shape:", dim(data$feature_metadata), "{newline}")
    }}

    # Load embeddings
    embedding_files <- list.files(data_dir, pattern = "^embedding_.*{ext}$", full.names = TRUE)
    if (length(embedding_files) > 0) {{
        cat("Loading embeddings...{newline}")
        data$embeddings <- list()
        for (emb_file in embedding_files) {{
            emb_name <- gsub("^embedding_(.*){ext}$", "\\\\1", basename(emb_file))
            emb_data <- as.data.frame(fread(emb_file))
            rownames(emb_data) <- emb_data[[1]]
            emb_data[[1]] <- NULL
            data$embeddings[[emb_name]] <- emb_data
            cat("Loaded", emb_name, "embedding:", dim(emb_data), "{newline}")
        }}
    }}

    return(data)
}}

# Function to create Seurat object
create_seurat_object <- function(data, project = "CSV_Import") {{
    if (is.null(data$count_matrix)) {{
        stop("Count matrix is required")
    }}

    cat("Creating Seurat object...{newline}")

    # Transpose count matrix (Seurat expects genes x cells)
    count_matrix <- t(as.matrix(data$count_matrix))

    # Create Seurat object
    seurat_obj <- CreateSeuratObject(
        counts = count_matrix,
        project = project
    )

    # Add cell metadata
    if (!is.null(data$cell_metadata)) {{
        common_cells <- intersect(colnames(seurat_obj), rownames(data$cell_metadata))
        if (length(common_cells) > 0) {{
            seurat_obj@meta.data <- cbind(
                seurat_obj@meta.data,
                data$cell_metadata[common_cells, , drop = FALSE]
            )
            cat("Added cell metadata for", length(common_cells), "cells{newline}")
        }}
    }}

    # Add feature metadata
    if (!is.null(data$feature_metadata)) {{
        common_features <- intersect(rownames(seurat_obj), rownames(data$feature_metadata))
        if (length(common_features) > 0) {{
            seurat_obj@assays$RNA@meta.features <- data$feature_metadata[common_features, , drop = FALSE]
            cat("Added feature metadata for", length(common_features), "features{newline}")
        }}
    }}

    # Add embeddings
    if (!is.null(data$embeddings)) {{
        for (emb_name in names(data$embeddings)) {{
            emb_data <- data$embeddings[[emb_name]]
            common_cells <- intersect(colnames(seurat_obj), rownames(emb_data))

            if (length(common_cells) > 0) {{
                # Determine key prefix
                if (emb_name == "pca") {{
                    key_prefix <- "PC_"
                }} else if (emb_name == "umap") {{
                    key_prefix <- "UMAP_"
                }} else if (emb_name == "tsne") {{
                    key_prefix <- "tSNE_"
                }} else {{
                    key_prefix <- paste0(toupper(substr(emb_name, 1, 1)),
                                       substr(emb_name, 2, nchar(emb_name)), "_")
                }}

                seurat_obj[[emb_name]] <- CreateDimReducObject(
                    embeddings = as.matrix(emb_data[common_cells, , drop = FALSE]),
                    key = key_prefix,
                    assay = DefaultAssay(seurat_obj)
                )
                cat("Added", emb_name, "embeddings{newline}")
            }}
        }}
    }}

    return(seurat_obj)
}}

# Usage examples:
# data <- load_csv_data()
# seurat_obj <- create_seurat_object(data)
'''

    # Save scripts
    python_file = os.path.join(output_dir, "load_csv_data.py")
    with open(python_file, 'w') as f:
        f.write(python_script_content)

    r_file = os.path.join(output_dir, "load_csv_data.R")
    with open(r_file, 'w') as f:
        f.write(r_script_content)

    logger.info(f"Loading scripts saved: {python_file}, {r_file}")

def generate_conversion_summary(adata, output_dir, compress):
    """
    Generate a summary report of the conversion

    Args:
        adata: AnnData object
        output_dir (str): Output directory path
        compress (bool): Whether files are compressed
    """
    logger.info("Generating conversion summary")

    # Calculate file sizes
    total_size = 0
    file_info = []

    for file_path in Path(output_dir).glob("*"):
        if file_path.is_file():
            size_mb = file_path.stat().st_size / (1024**2)
            total_size += size_mb
            file_info.append(f"- {file_path.name}: {size_mb:.2f} MB")

    # Calculate basic statistics
    if sparse.issparse(adata.X):
        sparsity = (1 - adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1])) * 100
        total_counts = adata.X.sum()
    else:
        sparsity = (np.sum(adata.X == 0) / adata.X.size) * 100
        total_counts = np.sum(adata.X)

    ext = ".csv.gz" if compress else ".csv"
    compression_note = "gzip compressed" if compress else "uncompressed"

    summary_content = f"""# H5AD to CSV Conversion Summary (v2.0)

## Conversion Details
- **Conversion Date**: {pd.Timestamp.now()}
- **Input Format**: H5AD (AnnData)
- **Output Format**: CSV ({compression_note})
- **Output Directory**: {os.path.basename(output_dir)}
- **Version**: 2.0.0 (fixed f-string syntax + gzip compression)

## Data Overview
- **Cells**: {adata.shape[0]:,}
- **Features**: {adata.shape[1]:,}
- **Total Counts**: {total_counts:,}
- **Sparsity**: {sparsity:.2f}%

## Metadata Information
- **Cell Metadata Columns**: {len(adata.obs.columns)}
- **Feature Metadata Columns**: {len(adata.var.columns)}

## Dimensionality Reductions
{chr(10).join([f"- {key.replace('X_', '').upper()}: {adata.obsm[key].shape[1]} dimensions" for key in adata.obsm.keys()]) if adata.obsm.keys() else "- None found"}

## Output Files
{chr(10).join(file_info)}

**Total Size**: {total_size:.2f} MB

## File Descriptions

### count_matrix{ext}
- Main expression/accessibility matrix
- Rows: cells, Columns: features
- Format: CSV with cell names as index, feature names as columns

### cell_metadata{ext}
- Cell annotations and quality metrics
- Rows: cells, Columns: metadata attributes
- Contains all information from AnnData .obs

### feature_metadata{ext}
- Feature annotations and statistics
- Rows: features, Columns: metadata attributes
- Contains all information from AnnData .var

### embedding_*{ext}
- Dimensionality reduction coordinates
- Separate file for each embedding (PCA, UMAP, t-SNE, etc.)
- Rows: cells, Columns: embedding dimensions

## Loading Instructions

### Python
```python
# Load data
exec(open('load_csv_data.py').read())
data = load_csv_data(".")

# Recreate AnnData object
adata = create_anndata_object(data)

# Basic analysis
import scanpy as sc
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata)
```

### R
```r
# Load data
source("load_csv_data.R")
data <- load_csv_data(".")

# Create Seurat object
seurat_obj <- create_seurat_object(data)

# Basic analysis
library(Seurat)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
DimPlot(seurat_obj, reduction = "pca")
```

## Compression Benefits
{f"- **Compression**: Enabled (gzip)" if compress else "- **Compression**: Disabled"}
{f"- **Space Savings**: Approximately 70-90% reduction in file size" if compress else "- **File Size**: Full uncompressed CSV format"}
{f"- **Compatibility**: Works with all CSV-reading tools" if compress else "- **Compatibility**: Standard CSV format"}

## Format Advantages
- **Universal Compatibility**: Works with Excel, R, Python, and most analysis tools
- **Human Readable**: Easy to inspect and manually edit if needed
- **No Dependencies**: No special libraries required for reading
- **Archival Format**: Long-term data preservation and sharing

## Performance Notes
- CSV format requires more storage than specialized formats (HDF5, MTX)
- Loading time depends on file size and compression
- Consider using compressed format for large datasets
- Memory usage scales with dataset size

## Usage Tips
- Use compressed format (csv.gz) to save disk space
- Load only required files for memory efficiency
- Consider chunking for very large datasets
- Validate data integrity after loading

## Version 2.0 Improvements
- Fixed f-string syntax error that prevented script execution
- Added gzip compression support for reduced file sizes
- Enhanced memory estimation and warnings
- Improved error handling and logging
- Better file organization and naming

## Troubleshooting
- **Large File Sizes**: Use compression option or consider MTX format
- **Memory Errors**: Reduce max_size_mb parameter or use chunking
- **Loading Issues**: Check file paths and compression settings
- **Missing Data**: Verify all expected files are present

## Notes
- CSV format creates dense matrices, requiring more space than sparse formats
- Best suited for smaller datasets or when universal compatibility is needed
- For large sparse datasets, consider MTX or HDF5 formats instead
- Compression is highly recommended for ATAC-seq data due to sparsity
"""

    summary_file = os.path.join(output_dir, "conversion_summary.md")
    with open(summary_file, 'w') as f:
        f.write(summary_content)

    logger.info(f"Summary saved: {summary_file}")

def main():
    """Main function with command line interface"""
    parser = argparse.ArgumentParser(
        description='Convert H5AD (AnnData) format to CSV format (Fixed Version v2.0)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s input.h5ad                           # Convert to CSV with compression
  %(prog)s input.h5ad -o csv_output             # Specify output directory
  %(prog)s input.h5ad --no-compress             # Save uncompressed CSV files
  %(prog)s input.h5ad --max-size-mb 500         # Limit file size warnings

Version 2.0 Features:
  - Fixed f-string syntax error
  - Added gzip compression support (default)
  - Enhanced memory estimation
  - Improved error handling
        '''
    )

    parser.add_argument('input',
                       help='Input H5AD file path')
    parser.add_argument('-o', '--output',
                       help='Output directory path (default: input filename + _csv)')
    parser.add_argument('--no-compress',
                       action='store_true',
                       help='Save uncompressed CSV files (default: compress with gzip)')
    parser.add_argument('--max-size-mb',
                       type=int,
                       default=1000,
                       help='Maximum uncompressed file size warning threshold in MB (default: 1000)')
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

    # Determine compression setting
    compress = not args.no_compress

    # Perform conversion
    success = h5ad_to_csv(
        args.input,
        args.output,
        compress=compress,
        max_size_mb=args.max_size_mb
    )

    if success:
        comp_status = "compressed" if compress else "uncompressed"
        print(f"\n✓ Conversion completed successfully!")
        print(f"CSV files ({comp_status}) saved in: {args.output or (os.path.splitext(args.input)[0] + '_csv')}")
        print(f"Load data with: python load_csv_data.py or source('load_csv_data.R')")
        print(f"Version: 2.0.0 (syntax fix + gzip compression)")
        sys.exit(0)
    else:
        print(f"\n✗ Conversion failed!")
        print(f"Check error messages above for troubleshooting.")
        sys.exit(1)

if __name__ == "__main__":
    main()