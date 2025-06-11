# 1. Use Pytoolkit_GeneSymbolUniform in comand line

## 1.1 Check parameters
```bash
python Pytoolkit_GeneSymbolUniform.py --help
```
***output:***

    usage: Pytoolkit_GeneSymbolUniform.py [-h] [--input_path INPUT_PATH] [--ref_table_path REF_TABLE_PATH] [--gene_list_path GENE_LIST_PATH] [--output_dir OUTPUT_DIR]
                                        [--output_prefix OUTPUT_PREFIX] [--print_report] [--average_alias]

    GeneSymbolUniform parameters

    optional arguments:
    -h, --help            show this help message and exit
    --input_path INPUT_PATH
                            the path of the input h5ad file
    --ref_table_path REF_TABLE_PATH
                            the path of the reference table.
    --gene_list_path GENE_LIST_PATH
                            the path of the total gene list table.
    --output_dir OUTPUT_DIR
                            the path to save the output h5ad file
    --output_prefix OUTPUT_PREFIX
                            the prefix of the output file and report file
    --print_report        print a report of the modified genes in the report.csv under the output_dir
    --average_alias       if average the counts of the genes mapped to the same aprroved symbol

## 1.2 Run Pytoolkit_GeneSymbolUniform
```bash
python Pytoolkit_GeneSymbolUniform.py --input_path /nfs/public/cell_gpt_data/bone_marrow.h5ad --print_report --average_alias
```

***output:***

    =========Task Started!=========
    =========Loading Reference Table=========
    Finished
    =========Processing Reference Table=========
    Finished
    =========Loading Query Data=========
    The shape of query data is: (85159, 33538) 
    Print out first 5 genes in query data, in case something wrong happens in data loading: 
    ['MIR1302-2HG' 'FAM138A' 'OR4F5' 'AL627309.1' 'AL627309.3']
    Finished
    =========Processing Gene List=========
    The length of reference gene_list is: 42117
    Finished
    =========Performing Gene Symbol Uniform=========
    Performing gene symbol uniform, this step may take several minutes
    Processing: 100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 33538/33538 [01:37<00:00, 342.92it/s]
    Finished
    =========Building output data=========
    Building output data, this step may take several minutes
    Processing: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 23360/23360 [00:22<00:00, 1045.87it/s]
    Shape of output data is (85159, 42117). It should have 42117 genes with cell number unchanged.
    =========Building output h5ad=========
    Finished
    =========Saving output h5ad=========
    h5ad file saved in: ./_uniformed.h5ad
    Finished
    =========Printing report=========
    report file saved in: ./_report.csv
    Finished


# 2. Use Pytoolkit_GeneSymbolUniform in Python

## 2.1 Check package
```python
from Pytoolkit_GeneSymbolUniform import GeneSymbolUniform
?GeneSymbolUniform
```

***output:***



    Signature:
    GeneSymbolUniform(
        input_adata=None,
        input_path=None,
        ref_table_path='./GeneSymbolRef_SelectAll_upd0731.csv',
        gene_list_path='./total_gene_list_42117.txt',
        output_dir='./',
        output_prefix='',
        print_report=True,
        average_alias=True,
    )

     Docstring:
        This function can be divided into following steps: 
        1) load reference table of the symbol relationships; 
        2) load the quert data (from a scanpy h5ad file with counts in 'X'); 
        3) construction the mapping dict between approved symbols and query symbols; 
        4) construct the output h5ad file and save the file and report
        

        :param input_adata: the input h5ad data with X as the count matrix to be uniformed.
        :type input_adata: scanpy::AnnData
        :default input_path: None

        :param input_path: the path of the input h5ad file (only used when input_data was not given)
        :type input_path: str
        :default input_path: None

        :param ref_table_path: the path of the reference table
        :type ref_table_path: str
        :default ref_table_path: "./GeneSymbolRef_SelectAll_upd0731.csv"

        :param gene_list_path: the path of the total gene list table
        :type gene_list_path: str
        :default gene_list_path: "./total_gene_list_42117.txt"

        :param output_dir: the path to save the output h5ad file
        :type output_dir: str
        :default output_dir: "./"

        :param output_prefix: the prefix of the output and report files
        :type output_prefix: str
        :default output_prefix: ""

        :param print_report: if print a report of the modified genes in the report.csv under the output_dir
        :type print_report: bool
        :defaul print_report: True

        :param average_alias: if average the counts of the genes mapped to the same aprroved symbol
        :type average_alias: bool
        :default average_alias: True


        :return: a h5ad data with uniformed epxression matrix.
        :rtype: scanpy::AnnData

     File:      ~/GeneSymbolUniform_Pytoolkit/Pytoolkit_GeneSymbolUniform.py
     Type:      function


## 2.2 Run Pytoolkit_GeneSymbolUniform

```python
import time
start_time = time.time()
adata_output = GeneSymbolUniform(input_path = "/nfs/public/cell_gpt_data/bone_marrow.h5ad",
                                 ref_table_path = "./GeneSymbolRef_SelectAll_upd0731.csv",
                                 gene_list_path = "total_gene_list_42117.txt",
                                 output_dir="./")
end_time = time.time()
```

***output:***


    =========Task Started!=========
    =========Loading Reference Table=========
    Finished
    =========Processing Reference Table=========
    Finished
    =========Loading Query Data=========
    The shape of query data is: (85159, 33538) 
    Print out first 5 genes in query data, in case something wrong happens in data loading: 
    ['MIR1302-2HG' 'FAM138A' 'OR4F5' 'AL627309.1' 'AL627309.3']
    Finished
    =========Processing Gene List=========
    The length of reference gene_list is: 42117
    Finished
    =========Performing Gene Symbol Uniform=========
    Performing gene symbol uniform, this step may take several minutes


    Processing: 100%|██████████| 33538/33538 [01:47<00:00, 312.78it/s]


    Finished
    =========Building output data=========
    Building output data, this step may take several minutes


    Processing: 100%|██████████| 23360/23360 [00:21<00:00, 1110.87it/s]


    Shape of output data is (85159, 42117). It should have 42117 genes with cell number unchanged.
    =========Building output h5ad=========
    Finished
    =========Saving output h5ad=========
    h5ad file saved in: ./_uniformed.h5ad
    Finished
    =========Printing report=========
    report file saved in: ./_report.csv
    Finished



```python
end_time-start_time
```

***output:***


    228.49409699440002


## 2.3 Check output


```python
adata_output
```


***output:***

    AnnData object with n_obs × n_vars = 85159 × 42117
        obs: 'original_name', 'cell_type', 'organ', 'sample_status', 'donor_age', 'age_bin', 'study_id', 'donor_gender'


