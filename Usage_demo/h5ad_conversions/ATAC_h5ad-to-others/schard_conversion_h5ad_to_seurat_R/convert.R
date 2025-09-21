#!/usr/bin/env Rscript

# H5AD to Seurat via schard Package (R Version v2.0)
#
# This script demonstrates conversion between H5AD format and Seurat objects
# using the schard R package, which provides direct H5AD file reading capabilities.
#
# Fixed: Correct function names from schard package
#
# Author: ATAC Analysis Team
# Version: 2.0.0

# Load required libraries
suppressPackageStartupMessages({
    library(optparse)
    library(tools)
})

# Function to print messages with timestamp
log_message <- function(message, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, message))
}

# Function to check if required packages are installed
check_requirements <- function() {
    log_message("Checking R package requirements")

    required_packages <- c("schard", "Seurat", "Signac")
    missing_packages <- c()

    for (package in required_packages) {
        if (!requireNamespace(package, quietly = TRUE)) {
            missing_packages <- c(missing_packages, package)
        }
    }

    if (length(missing_packages) > 0) {
        log_message(sprintf("Missing required R packages: %s", paste(missing_packages, collapse = ", ")), "ERROR")
        log_message("Install missing packages with:", "ERROR")
        for (pkg in missing_packages) {
            if (pkg == "schard") {
                log_message("  remotes::install_github('cellgeni/schard')", "ERROR")
            } else {
                log_message(sprintf("  install.packages('%s')", pkg), "ERROR")
            }
        }
        return(FALSE)
    }

    log_message("All required R packages are available")
    return(TRUE)
}

# Function to validate input file
validate_input_file <- function(h5ad_file) {
    if (!file.exists(h5ad_file)) {
        log_message(sprintf("Input file does not exist: %s", h5ad_file), "ERROR")
        return(FALSE)
    }

    if (!grepl("\\.h5ad$", h5ad_file)) {
        log_message(sprintf("File does not have .h5ad extension: %s", h5ad_file), "WARNING")
    }

    return(TRUE)
}

# Function to convert H5AD to Seurat using schard
convert_h5ad_to_seurat_schard <- function(h5ad_file, output_file = NULL, use_raw = FALSE, spatial = FALSE, simplify_spatial = TRUE) {
    log_message(sprintf("Converting %s to Seurat object using schard", h5ad_file))

    # Set default output file in current directory
    if (is.null(output_file)) {
        output_file <- file.path(getwd(), paste0(tools::file_path_sans_ext(basename(h5ad_file)), "_seurat_schard.rds"))
    }

    # Create output directory in current working directory
    output_dir <- getwd()
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        log_message(sprintf("Created output directory: %s", output_dir))
    }

    tryCatch({
        # Load required libraries
        library(schard)
        library(Seurat)
        library(Signac)

        # Convert based on data type
        log_message("Loading H5AD file using schard")

        if (spatial) {
            log_message("Converting spatial data with h5ad2seurat_spatial")
            seurat_obj <- schard::h5ad2seurat_spatial(
                h5ad_file,
                use.raw = use_raw,
                simplify = simplify_spatial
            )
        } else {
            log_message("Converting regular data with h5ad2seurat")
            seurat_obj <- schard::h5ad2seurat(
                h5ad_file,
                use.raw = use_raw
            )
        }

        # Basic information about the converted object
        if (spatial && !simplify_spatial) {
            # List of Seurat objects
            log_message(sprintf("Converted to list with %d Seurat objects", length(seurat_obj)))
            for (i in seq_along(seurat_obj)) {
                obj <- seurat_obj[[i]]
                log_message(sprintf("Object %s: %d cells, %d features",
                                  names(seurat_obj)[i], ncol(obj), nrow(obj)))
            }
        } else {
            # Single Seurat object
            log_message(sprintf("Object summary: %d cells, %d features", ncol(seurat_obj), nrow(seurat_obj)))

            # Print assay information
            assays <- names(seurat_obj@assays)
            log_message(sprintf("Available assays: %s", paste(assays, collapse = ", ")))

            # Print reduction information
            if (length(seurat_obj@reductions) > 0) {
                reductions <- names(seurat_obj@reductions)
                log_message(sprintf("Available reductions: %s", paste(reductions, collapse = ", ")))
            } else {
                log_message("No dimensionality reductions found")
            }

            # Print spatial information if available
            if (length(seurat_obj@images) > 0) {
                images <- names(seurat_obj@images)
                log_message(sprintf("Available spatial images: %s", paste(images, collapse = ", ")))
            }
        }

        log_message("Seurat object created successfully")

        # Save Seurat object
        log_message(sprintf("Saving Seurat object to: %s", output_file))
        saveRDS(seurat_obj, file = output_file)

        log_message("Conversion completed successfully!")
        return(TRUE)

    }, error = function(e) {
        log_message(sprintf("Conversion failed: %s", e$message), "ERROR")
        return(FALSE)
    })
}

# Function to convert H5AD to SingleCellExperiment using schard
convert_h5ad_to_sce_schard <- function(h5ad_file, output_file = NULL) {
    log_message(sprintf("Converting %s to SingleCellExperiment object using schard", h5ad_file))

    # Set default output file in current directory
    if (is.null(output_file)) {
        output_file <- file.path(getwd(), paste0(tools::file_path_sans_ext(basename(h5ad_file)), "_sce_schard.rds"))
    }

    # Create output directory in current working directory
    output_dir <- getwd()
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        log_message(sprintf("Created output directory: %s", output_dir))
    }

    tryCatch({
        # Load required libraries
        library(schard)
        library(SingleCellExperiment)

        # Convert using schard
        log_message("Loading H5AD file as SingleCellExperiment using schard")
        sce_obj <- schard::h5ad2sce(h5ad_file)

        # Basic information
        log_message(sprintf("SCE object summary: %d cells, %d features", ncol(sce_obj), nrow(sce_obj)))

        # Print assay information
        assays <- assayNames(sce_obj)
        log_message(sprintf("Available assays: %s", paste(assays, collapse = ", ")))

        # Print reduction information
        if (length(reducedDims(sce_obj)) > 0) {
            reductions <- names(reducedDims(sce_obj))
            log_message(sprintf("Available reduced dimensions: %s", paste(reductions, collapse = ", ")))
        } else {
            log_message("No dimensionality reductions found")
        }

        log_message("SingleCellExperiment object created successfully")

        # Save SCE object
        log_message(sprintf("Saving SingleCellExperiment object to: %s", output_file))
        saveRDS(sce_obj, file = output_file)

        log_message("Conversion completed successfully!")
        return(TRUE)

    }, error = function(e) {
        log_message(sprintf("Conversion failed: %s", e$message), "ERROR")
        return(FALSE)
    })
}

# Function to create usage examples
create_usage_examples <- function(output_dir) {
    log_message("Creating usage examples")

    # R usage example with schard workflows
    r_example_content <- '# R Usage Examples for schard-based H5AD Conversion (v2.0)
# Generated by schard conversion script v2.0

library(schard)
library(Seurat)
library(Signac)
library(SingleCellExperiment)

# Example 1: Basic H5AD to Seurat conversion
convert_basic_seurat <- function(h5ad_file) {
    cat("Converting H5AD to Seurat object using schard\\n")

    # Load as Seurat object
    seurat_obj <- schard::h5ad2seurat(h5ad_file)

    # Basic inspection
    cat("Object summary:\\n")
    print(seurat_obj)

    # Check available assays
    cat("Available assays:", paste(names(seurat_obj@assays), collapse = ", "), "\\n")

    # Check reductions
    if (length(seurat_obj@reductions) > 0) {
        cat("Available reductions:", paste(names(seurat_obj@reductions), collapse = ", "), "\\n")

        # Plot first available reduction
        reduction_name <- names(seurat_obj@reductions)[1]
        cat("Plotting", reduction_name, "\\n")
        print(DimPlot(seurat_obj, reduction = reduction_name))
    }

    return(seurat_obj)
}

# Example 2: Load raw counts instead of processed
convert_raw_seurat <- function(h5ad_file) {
    cat("Converting H5AD to Seurat with raw counts\\n")

    # Load raw counts
    seurat_raw <- schard::h5ad2seurat(h5ad_file, use.raw = TRUE)

    # Compare with processed counts
    seurat_processed <- schard::h5ad2seurat(h5ad_file, use.raw = FALSE)

    cat("Raw counts summary:\\n")
    print(seurat_raw)
    cat("Processed counts summary:\\n")
    print(seurat_processed)

    # Compare count totals
    raw_totals <- colSums(seurat_raw)
    proc_totals <- colSums(seurat_processed)

    cat("Raw vs processed count correlation:", cor(raw_totals, proc_totals), "\\n")

    return(list(raw = seurat_raw, processed = seurat_processed))
}

# Example 3: Spatial data conversion
convert_spatial_data <- function(h5ad_file) {
    cat("Converting spatial H5AD to Seurat\\n")

    # Option 1: Single merged object
    spatial_merged <- schard::h5ad2seurat_spatial(h5ad_file, simplify = TRUE)

    # Option 2: List of objects per slide
    spatial_list <- schard::h5ad2seurat_spatial(h5ad_file, simplify = FALSE)

    cat("Merged spatial object:\\n")
    print(spatial_merged)

    cat("List of spatial objects:\\n")
    print(sapply(spatial_list, function(x) paste(ncol(x), "cells")))

    # Spatial plotting
    if (length(spatial_merged@images) > 0) {
        cat("Creating spatial plots\\n")

        # Plot total counts
        spatial_plot <- Seurat::SpatialPlot(spatial_merged, features = "nCount_RNA")
        print(spatial_plot)

        # Plot for specific image if multiple
        image_names <- names(spatial_merged@images)
        if (length(image_names) > 1) {
            specific_plot <- Seurat::SpatialPlot(
                spatial_merged,
                features = "nCount_RNA",
                images = image_names[1]
            )
            print(specific_plot)
        }
    }

    return(list(merged = spatial_merged, list = spatial_list))
}

# Example 4: SingleCellExperiment conversion
convert_to_sce <- function(h5ad_file) {
    cat("Converting H5AD to SingleCellExperiment\\n")

    # Load as SCE
    sce_obj <- schard::h5ad2sce(h5ad_file)

    # Basic inspection
    cat("SCE object summary:\\n")
    print(sce_obj)

    # Check assays
    cat("Available assays:", paste(assayNames(sce_obj), collapse = ", "), "\\n")

    # Check reduced dimensions
    if (length(reducedDims(sce_obj)) > 0) {
        cat("Available reduced dimensions:", paste(names(reducedDims(sce_obj)), collapse = ", "), "\\n")

        # Plot first reduction
        rd_name <- names(reducedDims(sce_obj))[1]
        rd_data <- reducedDim(sce_obj, rd_name)

        if (ncol(rd_data) >= 2) {
            plot(rd_data[,1], rd_data[,2],
                 main = paste(rd_name, "plot"),
                 xlab = paste(rd_name, "1"),
                 ylab = paste(rd_name, "2"),
                 pch = 16, cex = 0.5)
        }
    }

    return(sce_obj)
}

# Example 5: ATAC-seq specific analysis
analyze_atac_data <- function(seurat_obj) {
    cat("Analyzing ATAC-seq data\\n")

    # Check if this looks like ATAC data
    feature_names <- rownames(seurat_obj)
    is_atac <- mean(grepl(":", feature_names)) > 0.5

    if (is_atac) {
        cat("Detected ATAC-seq data (genomic coordinates)\\n")

        # Basic ATAC QC if not already present
        if (!"nCount_ATAC" %in% colnames(seurat_obj@meta.data)) {
            seurat_obj$nCount_ATAC <- Matrix::colSums(seurat_obj@assays[[DefaultAssay(seurat_obj)]]@counts)
            seurat_obj$nFeature_ATAC <- Matrix::colSums(seurat_obj@assays[[DefaultAssay(seurat_obj)]]@counts > 0)
        }

        # Basic filtering
        cat("Applying basic ATAC-seq filtering\\n")
        seurat_obj <- subset(
            seurat_obj,
            subset = nCount_ATAC > 1000 &
                     nCount_ATAC < 50000 &
                     nFeature_ATAC > 500
        )

        cat("After filtering:", ncol(seurat_obj), "cells remain\\n")

        # Standard ATAC workflow if not already done
        if (!"lsi" %in% names(seurat_obj@reductions)) {
            cat("Running standard ATAC-seq workflow\\n")

            # Use Signac for ATAC processing
            library(Signac)

            # Rename assay to ATAC if needed
            if (DefaultAssay(seurat_obj) != "ATAC") {
                seurat_obj[["ATAC"]] <- seurat_obj[[DefaultAssay(seurat_obj)]]
                DefaultAssay(seurat_obj) <- "ATAC"
            }

            # TF-IDF normalization
            seurat_obj <- RunTFIDF(seurat_obj)

            # Find top features
            seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = "q0")

            # SVD/LSI
            seurat_obj <- RunSVD(seurat_obj)

            # Clustering
            seurat_obj <- FindNeighbors(seurat_obj, reduction = "lsi", dims = 2:30)
            seurat_obj <- FindClusters(seurat_obj, algorithm = 3)

            # UMAP
            seurat_obj <- RunUMAP(seurat_obj, reduction = "lsi", dims = 2:30)

            cat("Completed standard ATAC-seq workflow\\n")
        }
    } else {
        cat("Data does not appear to be ATAC-seq (no genomic coordinates detected)\\n")
    }

    return(seurat_obj)
}

# Example 6: Quality control and comparison
quality_control_analysis <- function(h5ad_file) {
    cat("Performing quality control analysis\\n")

    # Load both raw and processed
    seurat_raw <- schard::h5ad2seurat(h5ad_file, use.raw = TRUE)
    seurat_proc <- schard::h5ad2seurat(h5ad_file, use.raw = FALSE)

    # Compare basic statistics
    raw_counts <- Matrix::colSums(seurat_raw)
    proc_counts <- Matrix::colSums(seurat_proc)

    raw_features <- Matrix::colSums(seurat_raw > 0)
    proc_features <- Matrix::colSums(seurat_proc > 0)

    # Create comparison plots
    par(mfrow = c(2, 2))

    # Count comparison
    plot(raw_counts, proc_counts,
         main = "Raw vs Processed Counts",
         xlab = "Raw Counts", ylab = "Processed Counts",
         pch = 16, cex = 0.5)
    abline(0, 1, col = "red", lty = 2)

    # Feature comparison
    plot(raw_features, proc_features,
         main = "Raw vs Processed Features",
         xlab = "Raw Features", ylab = "Processed Features",
         pch = 16, cex = 0.5)
    abline(0, 1, col = "red", lty = 2)

    # Distribution plots
    hist(log10(raw_counts + 1), main = "Raw Counts Distribution",
         xlab = "log10(Counts + 1)", breaks = 50)
    hist(log10(proc_counts + 1), main = "Processed Counts Distribution",
         xlab = "log10(Counts + 1)", breaks = 50)

    par(mfrow = c(1, 1))

    # Summary statistics
    cat("\\nQuality Control Summary:\\n")
    cat("Raw data - Median counts per cell:", median(raw_counts), "\\n")
    cat("Processed data - Median counts per cell:", median(proc_counts), "\\n")
    cat("Raw data - Median features per cell:", median(raw_features), "\\n")
    cat("Processed data - Median features per cell:", median(proc_features), "\\n")
    cat("Count correlation:", cor(raw_counts, proc_counts), "\\n")

    return(list(
        raw = seurat_raw,
        processed = seurat_proc,
        raw_counts = raw_counts,
        proc_counts = proc_counts
    ))
}

# Usage examples workflow:
# 1. Basic conversion
# seurat_obj <- convert_basic_seurat("data.h5ad")

# 2. Raw vs processed comparison
# comparison <- convert_raw_seurat("data.h5ad")

# 3. Spatial data analysis
# spatial_results <- convert_spatial_data("spatial_data.h5ad")

# 4. SingleCellExperiment conversion
# sce_obj <- convert_to_sce("data.h5ad")

# 5. ATAC-seq specific analysis
# atac_obj <- analyze_atac_data(seurat_obj)

# 6. Quality control
# qc_results <- quality_control_analysis("data.h5ad")
'

    # Save example script
    example_file <- file.path(output_dir, "schard_usage_examples.R")
    writeLines(r_example_content, example_file)

    log_message(sprintf("Usage examples saved to: %s", example_file))
}

# Function to generate conversion summary
generate_conversion_summary <- function(h5ad_file, output_file, conversion_type) {
    log_message("Generating conversion summary")

    output_dir <- getwd()

    # Load the saved object to get info
    converted_obj <- readRDS(output_file)

    # Get object information based on type
    if (inherits(converted_obj, "Seurat")) {
        obj_type <- "Seurat"
        n_cells <- ncol(converted_obj)
        n_features <- nrow(converted_obj)
        assays <- paste(names(converted_obj@assays), collapse = ", ")
        reductions <- paste(names(converted_obj@reductions), collapse = ", ")
        spatial_info <- if (length(converted_obj@images) > 0) {
            paste("Spatial images:", paste(names(converted_obj@images), collapse = ", "))
        } else {
            "No spatial data"
        }
    } else if (inherits(converted_obj, "SingleCellExperiment")) {
        obj_type <- "SingleCellExperiment"
        n_cells <- ncol(converted_obj)
        n_features <- nrow(converted_obj)
        assays <- paste(assayNames(converted_obj), collapse = ", ")
        reductions <- paste(names(reducedDims(converted_obj)), collapse = ", ")
        spatial_info <- "N/A for SCE"
    } else if (is.list(converted_obj)) {
        obj_type <- "List of Seurat objects"
        n_cells <- sum(sapply(converted_obj, ncol))
        n_features <- "Variable"
        assays <- "Multiple objects"
        reductions <- "Multiple objects"
        spatial_info <- sprintf("Contains %d spatial slides", length(converted_obj))
    }

    summary_content <- sprintf("# H5AD to %s Conversion via schard (R Version v2.0)

## Conversion Details
- **Date**: %s
- **Method**: schard package v2.0 (direct H5AD reading)
- **Input**: %s
- **Output**: %s
- **Conversion Type**: %s

## Data Overview
- **Object Type**: %s
- **Cells**: %s
- **Features**: %s
- **Assays**: %s
- **Reductions**: %s
- **Spatial Info**: %s

## Version 2.0 Improvements

### Correct Function Names
- **h5ad2seurat()**: Convert H5AD to Seurat object
- **h5ad2seurat_spatial()**: Convert spatial H5AD to Seurat with spatial information
- **h5ad2sce()**: Convert H5AD to SingleCellExperiment object
- **Raw data support**: use.raw parameter for accessing raw counts
- **Spatial options**: simplify parameter for merged vs list output

### schard Package Features
- **Direct H5AD reading**: No Python dependencies required
- **Native R integration**: Full compatibility with Seurat and SCE workflows
- **Spatial data support**: Handles Visium and other spatial technologies
- **Raw/processed data**: Access to both raw and normalized counts
- **Automatic metadata transfer**: Cell and feature annotations preserved

## About schard v2.0

The schard package provides native R functions for reading H5AD files without requiring
Python or reticulate. Version 2.0 uses the correct function names from the package.

### Key Advantages
- **No Python required**: Pure R implementation
- **Direct file access**: Fast H5AD reading without conversion
- **Spatial support**: Native handling of spatial transcriptomics data
- **Multiple output formats**: Seurat objects, SCE objects, or lists
- **Metadata preservation**: All annotations and embeddings transferred

## Workflow Overview

1. **Install schard** from GitHub (cellgeni/schard)
2. **Load H5AD directly** using appropriate schard function
3. **Choose output format** (Seurat, SCE, spatial options)
4. **Access raw or processed** data as needed
5. **Continue with standard analysis** pipelines

## Usage Examples

### Basic Conversion
```r
library(schard)

# Load as Seurat
seurat_obj <- schard::h5ad2seurat('%s')

# Load raw counts
seurat_raw <- schard::h5ad2seurat('%s', use.raw = TRUE)

# Load as SingleCellExperiment
sce_obj <- schard::h5ad2sce('%s')
```

### Spatial Data
```r
# Single merged spatial object
spatial_merged <- schard::h5ad2seurat_spatial('spatial.h5ad')

# List of objects per slide
spatial_list <- schard::h5ad2seurat_spatial('spatial.h5ad', simplify = FALSE)

# Spatial plotting
Seurat::SpatialPlot(spatial_merged, features = 'total_counts')

# Plot specific slide
Seurat::SpatialPlot(spatial_list$slide1, features = 'total_counts')
```

### Analysis Examples
```r
# Basic inspection
print(seurat_obj)
str(seurat_obj, max.level = 2)

# Check reductions (note: scanpy reductions are prefixed, e.g., 'Xumap_')
DimPlot(seurat_obj, reduction = 'Xumap_')  # if UMAP from scanpy

# Compare raw vs processed
plot(colSums(seurat_raw), colSums(seurat_obj), pch=16)
```

## Requirements

### R Packages
```r
# Install schard from GitHub
remotes::install_github('cellgeni/schard')

# Install dependencies
install.packages(c('Seurat', 'Signac', 'SingleCellExperiment'))
```

### System Requirements
- R >= 4.0.0
- No Python dependencies
- HDF5 libraries (usually included with R packages)

## Troubleshooting v2.0

### Function Name Issues (FIXED)
- **v1.0 error**: `'h5ad_to_seurat' is not an exported object`
- **v2.0 solution**: Use correct function names: `h5ad2seurat`, `h5ad2seurat_spatial`, `h5ad2sce`

### Common Issues
- **Package not found**: Install with `remotes::install_github('cellgeni/schard')`
- **File not found**: Check H5AD file path and permissions
- **Memory issues**: For large files, consider using raw=TRUE or processing subsets

### Reduction Names
- scanpy reductions are auto-translated with prefixes (e.g., 'Xumap_', 'Xpca_')
- Specify reduction explicitly in plotting: `DimPlot(obj, reduction = 'Xumap_')`
- Use `names(seurat_obj@reductions)` to see available reductions

## Performance Notes
- **Fast loading**: Direct HDF5 access without Python overhead
- **Memory efficient**: No intermediate format conversion
- **Spatial support**: Efficient handling of spatial coordinates and images
- **Large files**: Better performance than Python-based approaches for large datasets

## Advantages of schard v2.0
- **Pure R solution**: No Python dependencies or environment setup
- **Correct function names**: Uses official schard API
- **Spatial support**: Native handling of spatial transcriptomics
- **Multiple formats**: Seurat and SingleCellExperiment output options
- **Raw data access**: Easy access to both raw and processed counts
- **Fast performance**: Direct HDF5 reading optimized for R

## Notes
- schard v2.0 provides the most straightforward H5AD to R conversion
- No Python environment setup required
- Spatial data support makes it ideal for Visium and similar technologies
- All metadata, embeddings, and spatial information preserved
- Compatible with standard Seurat and SingleCellExperiment workflows
",
    obj_type,
    format(Sys.time(), "%%Y-%%m-%%d %%H:%%M:%%S"),
    basename(h5ad_file),
    basename(output_file),
    conversion_type,
    obj_type,
    format(n_cells, big.mark = ","),
    format(n_features, big.mark = ","),
    assays,
    reductions,
    spatial_info,
    basename(h5ad_file),
    basename(h5ad_file),
    basename(h5ad_file)
)

    summary_file <- file.path(output_dir, "conversion_summary.md")
    writeLines(summary_content, summary_file)
    log_message(sprintf("Summary saved to: %s", summary_file))
}

# Main function
main <- function() {
    # Define command line options
    option_list <- list(
        make_option(c("-i", "--input"), type = "character", default = NULL,
                   help = "Input H5AD file path", metavar = "character"),
        make_option(c("-o", "--output"), type = "character", default = NULL,
                   help = "Output RDS file path (default: input filename with _seurat_schard.rds suffix)",
                   metavar = "character"),
        make_option(c("--use-raw"), action = "store_true", default = FALSE,
                   help = "Use raw counts instead of processed data"),
        make_option(c("--spatial"), action = "store_true", default = FALSE,
                   help = "Convert spatial data using h5ad2seurat_spatial"),
        make_option(c("--no-simplify"), action = "store_true", default = FALSE,
                   help = "For spatial data, return list of objects instead of merged object"),
        make_option(c("--to-sce"), action = "store_true", default = FALSE,
                   help = "Convert to SingleCellExperiment instead of Seurat"),
        make_option(c("--check-only"), action = "store_true", default = FALSE,
                   help = "Only check if requirements are installed"),
        make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                   help = "Enable verbose logging")
    )

    # Parse command line arguments
    opt_parser <- OptionParser(option_list = option_list,
                              description = "Convert H5AD to Seurat/SCE object using schard package v2.0",
                              epilogue = "Examples:
  Rscript convert.R -i input.h5ad                           # Convert to Seurat
  Rscript convert.R -i input.h5ad --use-raw                 # Use raw counts
  Rscript convert.R -i input.h5ad --spatial                 # Spatial data (merged)
  Rscript convert.R -i input.h5ad --spatial --no-simplify   # Spatial data (list)
  Rscript convert.R -i input.h5ad --to-sce                  # Convert to SCE
  Rscript convert.R --check-only                            # Check requirements")

    opt <- parse_args(opt_parser)

    # Check R requirements
    if (!check_requirements()) {
        quit(status = 1)
    }

    if (opt$`check-only`) {
        cat("\n✓ All requirements are satisfied!\n")
        cat("You can now convert H5AD files using schard v2.0.\n")
        cat("\nAvailable functions:\n")
        cat("- schard::h5ad2seurat() - Convert to Seurat object\n")
        cat("- schard::h5ad2seurat_spatial() - Convert spatial data\n")
        cat("- schard::h5ad2sce() - Convert to SingleCellExperiment\n")
        quit(status = 0)
    }

    # Check if input file is provided
    if (is.null(opt$input)) {
        print_help(opt_parser)
        quit(status = 1)
    }

    # Validate input file
    if (!validate_input_file(opt$input)) {
        quit(status = 1)
    }

    # Set output file in current directory
    conversion_type <- "Seurat"
    if (opt$`to-sce`) {
        conversion_type <- "SingleCellExperiment"
        if (is.null(opt$output)) {
            opt$output <- file.path(getwd(), paste0(tools::file_path_sans_ext(basename(opt$input)), "_sce_schard.rds"))
        }
    } else if (opt$spatial) {
        conversion_type <- if (opt$`no-simplify`) "Spatial Seurat List" else "Spatial Seurat"
        if (is.null(opt$output)) {
            opt$output <- file.path(getwd(), paste0(tools::file_path_sans_ext(basename(opt$input)), "_spatial_schard.rds"))
        }
    } else {
        if (is.null(opt$output)) {
            opt$output <- file.path(getwd(), paste0(tools::file_path_sans_ext(basename(opt$input)), "_seurat_schard.rds"))
        }
    }

    # Create output directory and examples in current directory
    create_usage_examples(getwd())

    # Perform conversion based on type
    if (opt$`to-sce`) {
        success <- convert_h5ad_to_sce_schard(opt$input, opt$output)
    } else {
        success <- convert_h5ad_to_seurat_schard(
            opt$input,
            opt$output,
            use_raw = opt$`use-raw`,
            spatial = opt$spatial,
            simplify_spatial = !opt$`no-simplify`
        )
    }

    if (success) {
        # Generate summary
        tryCatch({
            generate_conversion_summary(opt$input, opt$output, conversion_type)
        }, error = function(e) {
            log_message("Could not generate summary", "WARNING")
        })

        cat(sprintf("\n✓ Conversion completed successfully!\n"))
        cat(sprintf("%s object saved to: %s\n", conversion_type, opt$output))
        cat(sprintf("Usage examples created in: %s\n", getwd()))
        cat(sprintf("Load in R with: obj <- readRDS('%s')\n", opt$output))
        cat("\nThe converted object is ready for analysis with Seurat/Signac workflows.\n")
        quit(status = 0)
    } else {
        cat(sprintf("\n✗ Conversion failed!\n"))
        cat("Check the error messages above for troubleshooting.\n")
        quit(status = 1)
    }
}

# Run main function if script is executed directly
if (!interactive()) {
    main()
}