# hECA-ATAC dataTrans Conversion Guide

## Python Conversion Scripts

### 1. h5ad_to_mtx_py# Create and activate environment
conda create -n heca-atac
conda activate heca-atac

# Navigate to directory and install dependencies
cd h5ad_to_mtx_py
pip install -r requirements.txt 

# Run conversion
python convert.py ../ATAC-10.1126%2Fscience.aba7612-Spleen.h5ad 
### 2. h5ad_to_loom_py_v4# Create and activate environment
conda create -n heca-atac
conda activate heca-atac

# Navigate to directory and install dependencies
cd h5ad_to_loom_py
pip install -r requirements.txt 

# Run conversion
python convert.py ../ATAC-10.1126%2Fscience.aba7612-Spleen.h5ad 
### 3. h5ad_to_seurat_py# Create and activate environment
conda create -n heca-atac
conda activate heca-atac

# Navigate to directory and install dependencies
cd h5ad_to_seurat_py
pip install -r requirements.txt 

# Run conversion
python convert.py ../ATAC-10.1126%2Fscience.aba7612-Spleen.h5ad 
### 4. h5ad_to_csv_py# Create and activate environment
conda create -n heca-atac
conda activate heca-atac

# Navigate to directory and install dependencies
cd h5ad_to_csv_py
pip install -r requirements.txt 

# (Conversion command not specified)
## R Conversion Script

### schard_conversion_h5ad_to_seurat_R# Install required package
devtools::install_github('cellgeni/schard')

# Navigate to directory and run conversion
cd schard_conversion_R
Rscript convert.R -i ../ATAC-10.1126%2Fscience.aba7612-Spleen.h5ad