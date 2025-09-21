# H5AD to Seurat/SCE via schard Package (R Version v2.0)

This tool provides direct conversion between H5AD format and Seurat/SingleCellExperiment objects using the schard R package, which offers native H5AD file reading without Python dependencies.

## 🎯 Version 2.0 - Correct Function Names

This version fixes the function name errors and uses the correct schard API:

### Fixed Issues in v2.0
- **✅ Correct Functions**: Now uses `h5ad2seurat()`, `h5ad2seurat_spatial()`, `h5ad2sce()`
- **✅ Spatial Data Support**: Proper handling of spatial transcriptomics data
- **✅ Raw Data Access**: Support for both raw and processed counts
- **✅ Multiple Output Formats**: Seurat objects, SCE objects, or lists of objects

### Key Function Names (Fixed)
```r
# v1.0 Problem: Wrong function names
schard::h5ad_to_seurat()  # ❌ Does not exist

# v2.0 Solution: Correct function names
schard::h5ad2seurat()           # ✅ Convert to Seurat
schard::h5ad2seurat_spatial()   # ✅ Convert spatial data
schard::h5ad2sce()             # ✅ Convert to SCE
```

## Installation and Usage

### Quick Start

```bash
# Start R
R

# Install schard from GitHub
remotes::install_github('cellgeni/schard')

# Install dependencies
install.packages(c("Seurat", "Signac", "SingleCellExperiment"))

# Exit R and test conversion
cd schard_conversion_R_v2
Rscript convert.R ../ATAC-10.1126%2Fscience.aba7612-Spleen.h5ad
```

### Expected Success Output
```
[2025-09-20 XX:XX:XX] INFO: Loading H5AD file using schard
[2025-09-20 XX:XX:XX] INFO: Converting regular data with h5ad2seurat
[2025-09-20 XX:XX:XX] INFO: Object summary: XXXX cells, XXXX features
[2025-09-20 XX:XX:XX] INFO: Available assays: RNA, ATAC
[2025-09-20 XX:XX:XX] INFO: Available reductions: Xumap_, Xpca_
[2025-09-20 XX:XX:XX] INFO: Seurat object created successfully

✓ Conversion completed successfully!
```

## Features

### Native R Implementation
- **No Python required**: Pure R solution using schard package
- **Direct H5AD reading**: Fast access to H5AD files without conversion
- **Memory efficient**: Optimized HDF5 reading for large datasets
- **Full compatibility**: Works with all Seurat and SingleCellExperiment workflows

### Spatial Data Support
- **Merged objects**: Single Seurat object with all spatial slides
- **List output**: Separate Seurat objects per slide for detailed analysis
- **Spatial coordinates**: Automatic handling of spatial metadata
- **Image data**: Preservation of spatial images and coordinates

### Multiple Output Formats
- **Seurat objects**: Standard single-cell analysis format
- **SingleCellExperiment**: Bioconductor-compatible format
- **Raw/processed data**: Access to both raw and normalized counts
- **Spatial variants**: Specialized handling for spatial transcriptomics

## Command Line Interface

```bash
# Basic conversion to Seurat
Rscript convert.R -i input.h5ad

# Use raw counts instead of processed
Rscript convert.R -i input.h5ad --use-raw

# Convert to SingleCellExperiment
Rscript convert.R -i input.h5ad --to-sce

# Spatial data conversion (merged object)
Rscript convert.R -i spatial.h5ad --spatial

# Spatial data conversion (list of objects)
Rscript convert.R -i spatial.h5ad --spatial --no-simplify

# Specify output file
Rscript convert.R -i input.h5ad -o my_seurat.rds

# Check requirements only
Rscript convert.R --check-only
```

## Usage Examples

### Basic Conversion
```r
library(schard)
library(Seurat)

# Load H5AD as Seurat object
seurat_obj <- schard::h5ad2seurat("data.h5ad")

# Check the object
print(seurat_obj)
DimPlot(seurat_obj, reduction = "Xumap_")  # Note: scanpy UMAP is "Xumap_"
```

### Raw vs Processed Data
```r
# Load processed data (default)
seurat_processed <- schard::h5ad2seurat("data.h5ad")

# Load raw counts
seurat_raw <- schard::h5ad2seurat("data.h5ad", use.raw = TRUE)

# Compare count totals
plot(colSums(seurat_raw), colSums(seurat_processed), pch=16)
```

### Spatial Data Analysis
```r
# Single merged spatial object
spatial_merged <- schard::h5ad2seurat_spatial("spatial.h5ad")

# List of objects per slide
spatial_list <- schard::h5ad2seurat_spatial("spatial.h5ad", simplify = FALSE)

# Spatial plotting
Seurat::SpatialPlot(spatial_merged, features = "total_counts")

# Plot specific slide
Seurat::SpatialPlot(spatial_merged, features = "total_counts",
                   images = "HCAHeartST11702009")

# Plot from list
Seurat::SpatialPlot(spatial_list$HCAHeartST11702010, features = "total_counts")
```

### SingleCellExperiment Conversion
```r
library(SingleCellExperiment)

# Load as SCE object
sce_obj <- schard::h5ad2sce("data.h5ad")

# Check the object
print(sce_obj)
plotReducedDim(sce_obj, dimred = "Xumap_")  # If available
```

## Data Structure Handling

### Reduction Names
schard automatically translates scanpy reductions with prefixes:
- `X_umap` → `Xumap_`
- `X_pca` → `Xpca_`
- `X_tsne` → `Xtsne_`

Use these names in Seurat functions:
```r
# List available reductions
names(seurat_obj@reductions)

# Use correct reduction name
DimPlot(seurat_obj, reduction = "Xumap_")
```

### Metadata Preservation
- **Cell metadata**: All `.obs` columns preserved in `@meta.data`
- **Feature metadata**: All `.var` columns preserved in assay metadata
- **Embeddings**: All `.obsm` matrices converted to Seurat reductions
- **Spatial data**: Coordinates and images properly handled

## Performance and Memory

### Advantages of schard
- **Fast loading**: Direct HDF5 access, no Python overhead
- **Memory efficient**: Optimized for large single-cell datasets
- **No dependencies**: Pure R solution, no Python environment setup
- **Spatial optimized**: Efficient handling of spatial coordinates and images

### Performance Characteristics
- **Small datasets** (< 10K cells): Near-instantaneous loading
- **Medium datasets** (10K-100K cells): Fast loading, low memory overhead
- **Large datasets** (> 100K cells): Excellent performance, efficient memory usage
- **Spatial data**: Optimized for Visium and similar technologies

## Troubleshooting

### Fixed Issues ✅

**v1.0 Error (FIXED)**:
```
'h5ad_to_seurat' is not an exported object from 'namespace:schard'
```

**v2.0 Result**:
```r
# Correct function names
schard::h5ad2seurat()
schard::h5ad2seurat_spatial()
schard::h5ad2sce()
```

### Common Issues

**Package installation**:
```r
# Install from GitHub
remotes::install_github('cellgeni/schard')

# Check installation
library(schard)
```

**Reduction names**:
```r
# Check available reductions
names(seurat_obj@reductions)

# Use correct prefix for scanpy reductions
DimPlot(seurat_obj, reduction = "Xumap_")
```

**Large files**:
```r
# Use raw data for smaller memory footprint
seurat_raw <- schard::h5ad2seurat("large_file.h5ad", use.raw = TRUE)
```

## Advanced Usage

### ATAC-seq Data
```r
# Load ATAC data
atac_obj <- schard::h5ad2seurat("atac_data.h5ad")

# Standard ATAC workflow with Signac
library(Signac)
atac_obj <- RunTFIDF(atac_obj)
atac_obj <- FindTopFeatures(atac_obj, min.cutoff = 'q0')
atac_obj <- RunSVD(atac_obj)
atac_obj <- FindNeighbors(atac_obj, reduction = 'lsi', dims = 2:30)
atac_obj <- FindClusters(atac_obj, algorithm = 3)
atac_obj <- RunUMAP(atac_obj, reduction = 'lsi', dims = 2:30)
```

### Quality Control Analysis
```r
# Compare raw vs processed
raw_obj <- schard::h5ad2seurat("data.h5ad", use.raw = TRUE)
proc_obj <- schard::h5ad2seurat("data.h5ad", use.raw = FALSE)

# QC metrics
VlnPlot(proc_obj, features = c("nFeature_RNA", "nCount_RNA"))
FeatureScatter(proc_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```

### Batch Processing
```r
# Process multiple files
h5ad_files <- list.files(pattern = "\\.h5ad$")
seurat_objects <- lapply(h5ad_files, function(f) {
    schard::h5ad2seurat(f)
})
names(seurat_objects) <- gsub("\\.h5ad$", "", h5ad_files)
```

## Requirements

### R Packages
```r
# Core packages
remotes::install_github('cellgeni/schard')  # Main package
install.packages(c('Seurat', 'Signac'))     # Analysis packages

# Optional for SingleCellExperiment
BiocManager::install('SingleCellExperiment')
```

### System Requirements
- **R** >= 4.0.0 (recommended: >= 4.1.0)
- **No Python** dependencies required
- **HDF5 libraries** (usually included with R packages)
- **Memory**: Varies by dataset size, but generally efficient

## Format Advantages

### Advantages of schard Conversion
- ✅ **Pure R solution**: No Python environment setup
- ✅ **Fast performance**: Direct HDF5 reading optimized for R
- ✅ **Spatial support**: Native handling of spatial transcriptomics
- ✅ **Multiple formats**: Seurat and SingleCellExperiment output
- ✅ **Raw data access**: Easy switching between raw and processed
- ✅ **Full compatibility**: Works with all Seurat/Signac workflows

### When to Use schard
- ✅ Want pure R workflow without Python dependencies
- ✅ Need fast, reliable H5AD to Seurat conversion
- ✅ Working with spatial transcriptomics data
- ✅ Require both raw and processed data access
- ✅ Need SingleCellExperiment compatibility
- ✅ Want minimal setup and maximum reliability

## Notes

- **v2.0 uses correct function names** from the official schard package
- **No Python dependencies** required for conversion
- **Spatial data is fully supported** with proper coordinate and image handling
- **All metadata and embeddings are preserved** during conversion
- **scanpy reduction names are auto-translated** with appropriate prefixes
- **Both raw and processed data** can be accessed easily
- **Memory usage is optimized** for large single-cell datasets
- **Compatible with all standard Seurat/Signac workflows**