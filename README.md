# HIPC Dashboard Pipeline
## Overview
The HIPC Dashboard provides a web interface to the immune signatures curated as part of the HIPC Signatures II project (NIAID).
This respository, hipc-dashboard-pipeline, provides code and data to generate submission files for the HIPC Dashboard (http://hipc-dashboard.org/).

## Curated data
Curated data and various mapping and translation files are placed in ./source_data.
"HIPC Dashboard.xlsx" contains the primary curated data, with one tab for gene data and another for cell-type data.

## R scripts
R scripts are in ./src
"generate_HIPC_submissions.R" is the main script.  It expects to be called from its source location.
It has a number of optional settings at the beginning.  In particular, one can choose to run the script for the gene data or the cell-type data in the input Excel file, 

## Output files
### Submission files
If it does not already exist, a directory called "../submissions" will be created relative to ./src.
The submission files are tab-delimited.
Gene submission files will be placed in a subdirectory called "../submissions/hipc_gene/".
Cell-type submissions will be placed in a subdirectory called "../submissions/hipc_ctf/".

### Inspection files
As a convenience for easy data inspection, CSV versions of the files are also created and placed in subdirectories "../submissions/hipc_gene_csv" and "../submissions/hipc_ctf_csv". 

### Log files
An additional directory, "../logfiles" is created relative to the ./src directory.  After the script has run, this directory contains a number of files tracing the data transformations and final summary data, as well as "recreated_template" files that represent the "HIPC Dashboard.xlsx" sheets after all transformations.


