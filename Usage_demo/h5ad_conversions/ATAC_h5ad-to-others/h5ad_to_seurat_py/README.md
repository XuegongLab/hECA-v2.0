# H5AD to Seurat Conversion Tool

This tool converts ATAC-seq data from H5AD format (AnnData) to Seurat-compatible format, creating MTX matrix files and metadata files for seamless import into R/Seurat.

## Overview

The converter transforms single-cell ATAC-seq data stored in Python's AnnData format (.h5ad files) into a format that can be easily imported into R's Seurat ecosystem. This is particularly useful for researchers who want to leverage both Python's scanpy and R's Seurat/Signac packages in their analysis workflow.

## Features

- **Complete data preservation**: Converts count matrices, cell metadata, feature metadata, and dimensionality reductions
- **Seurat compatibility**: Generates files in the exact format expected by Seurat's data import functions
- **Automatic R script generation**: Creates a ready-to-use R script for importing the converted data
- **Data validation**: Includes comprehensive validation of ATAC-seq data format and quality
- **Memory efficient**: Preserves sparse matrix format to minimize memory usage
- **Detailed logging**: Provides comprehensive progress reporting and error handling

## Output Files

The conversion process generates the following files:

- `matrix.mtx`: Sparse count matrix in MTX format (features × cells, transposed for Seurat)
- `features.tsv`: Feature names and types (peak coordinates for ATAC data)
- `barcodes.tsv`: Cell barcodes/identifiers
- `cell_metadata.csv`: Cell annotations and metadata (if available)
- `feature_metadata.csv`: Feature annotations and metadata (if available)
- `*_embeddings.csv`: Dimensionality reduction results (PCA, UMAP, t-SNE, etc.)
- `import_seurat.R`: R script for importing the converted data into Seurat
- `conversion_summary.md`: Detailed summary of the conversion process

## Usage

### Command Line Interface

```bash
# Basic conversion
python convert.py input.h5ad

# Specify output directory
python convert.py input.h5ad -o seurat_data

# Enable verbose logging
python convert.py input.h5ad --verbose
```

### In R (after conversion)

```r
# Source the generated import script
source("import_seurat.R")

# Import the converted data
atac_data <- import_converted_data()

# Verify the import
print(atac_data)

# Visualize if embeddings are available
if ("umap" %in% names(atac_data@reductions)) {
    DimPlot(atac_data, reduction = "umap")
}
```

## Data Requirements

### Input Format
- **File type**: H5AD (AnnData format)
- **Data type**: ATAC-seq count data (preferably sparse)
- **Feature names**: Genomic coordinates in format `chr:start-end` (recommended)
- **Count data**: Non-negative integer counts

### Expected Structure
```
AnnData object with:
- .X: Count matrix (cells × features)
- .obs: Cell metadata
- .var: Feature metadata
- .obsm: Dimensionality reductions (optional)
```

## Installation Requirements

See `requirements.txt` for Python dependencies. The tool requires:

- Python 3.7+
- scanpy for H5AD file reading
- pandas for data manipulation
- scipy for sparse matrix operations
- numpy for numerical operations

For R usage after conversion:
- R 4.0+
- Seurat package
- Matrix package

## Validation and Quality Control

The tool performs several validation checks:

1. **File existence and format validation**
2. **Feature name format checking** (warns if not in genomic coordinate format)
3. **Data range validation** (checks for negative values)
4. **Sparsity analysis**
5. **Matrix format verification**

## Troubleshooting

### Common Issues

**File not found errors**:
- Verify the input H5AD file path is correct
- Ensure you have read permissions for the input file

**Memory errors**:
- The tool automatically converts to sparse format for memory efficiency
- For very large datasets, ensure sufficient RAM is available

**Import errors in R**:
- Verify Seurat and Matrix packages are installed
- Check that all output files are present in the specified directory

**Feature name warnings**:
- ATAC-seq features should be in genomic coordinate format (chr:start-end)
- This is a warning only and won't prevent conversion

### Performance Notes

- Conversion time scales with dataset size
- Sparse matrices are preserved for memory efficiency
- Large embeddings may take additional time to export

## Example Workflow

```bash
# Convert ATAC-seq data
python convert.py atac_data.h5ad -o converted_data --verbose

# In R session
source("converted_data/import_seurat.R")
atac_obj <- import_converted_data("converted_data")

# Basic Seurat workflow
library(Signac)
library(Seurat)

# Add gene annotations (optional)
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# Annotation(atac_obj) <- annotations

# Basic quality control and analysis
atac_obj <- RunTFIDF(atac_obj)
atac_obj <- FindTopFeatures(atac_obj, min.cutoff = 'q0')
atac_obj <- RunSVD(atac_obj)

# Visualization
DimPlot(atac_obj, reduction = "umap")  # if UMAP was included
```

## Notes

- The count matrix is transposed during conversion (features × cells) to match Seurat's expected format
- All metadata and dimensionality reductions from the original H5AD file are preserved
- The tool is optimized for ATAC-seq data but can handle other single-cell data types
- Sparse matrix format is maintained throughout the conversion process for memory efficiency