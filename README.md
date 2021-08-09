# HIPC Dashboard Pipeline
## Overview
The HIPC Dashboard provides a web interface to the immune signatures curated as part of the HIPC Signatures II project (NIAID).
This respository, hipc-dashboard-pipeline, provides code and data to generate submission files for the HIPC Dashboard (http://hipc-dashboard.org/).

In general, curators initially enter data exactly as it appears in a publication.  This data is then standardized as needed using suitable ontologies.  The original annotations are preserved for qualtity control and provenance.

Standardization varies according to data type:
* Update gene symbols.  Curated gene symbols are updated to current HGNC/NCBI symbols based on (1) NCBI synonyms and (2) a manually created mapping table.  The later deals with specific problematic symbols found in the curated data, where examination of the original data is able to support a definite mapping.
* Standardize cell types.  Cell types and their markers are initially curated as they appear in each publication.  We create a mapping table to transform the orignal values to standardized values using the Cell Ontology and Protein Ontology.

## Curated data
Curated data and various mapping and translation files are placed in ./source_data.
The primary input data processed by the R script pipeline is stored in tab-delimited text files, one for each type of curated data.  These are currently:
* HIPC Dashboard - Gene Expression.tsv
* HIPC Dashboard - Cell type Frequency.tsv

The same data is also available in an Excel file, "HIPC Dashboard.xlsx", for easy review, however this file is not used by the script.

## R scripts
R scripts are in ./src.
"generate_HIPC_submissions.R" is the main script.  It expects to be called from its source location.
It has a number of optional settings at the beginning.  In particular, one can choose to run the script for the gene data or the cell-type data.

## Output files
### Submission files
If it does not already exist, a directory called "../submissions" will be created relative to ./src.
The submission files are tab-delimited.
Gene submission files are placed in a subdirectory called "../submissions/hipc_gene/".
Cell-type submissions are placed in a subdirectory called "../submissions/hipc_ctf/".

For easy inspection, CSV-formatted versions of the same files are also generated under "../submissions/hipc_gene_csv/" and "../submissions/hipc_ctf_csv/"

### Log files
An additional directory (not checked-in to GitHub), "../logfiles" is created relative to the ./src directory.  After the script has run, this directory contains a number of files tracing the data transformations and final summary data, as well as "recreated_template" files that represent the data after all updates and transformations, in the same format as the original data.

### Reformatted data files
Certain files from the "logfiles" directory are checked in to GitHub in the "reformatted_data" directory.  These are the "recreated_template" files, with identical data presented in tab-delimited, RDS and Excel formats, as well as the response components for each signature in Broad GMT format (tab-delimited).  These files represent the data after all updates and transformations. The "recreated_template" files are in the same spreadsheet format as the original data.


## Database Schema
![Dashboard DB schema](hipc_dashboard_data_model_simple.png)
