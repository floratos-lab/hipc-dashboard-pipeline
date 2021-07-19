# HIPC Dashboard Pipeline
## Overview
The HIPC Dashboard provides a web interface to the immune signatures curated as part of the HIPC Signatures II project (NIAID).
This respository, hipc-dashboard-pipeline, provides code and data to generate submission files for the HIPC Dashboard (http://hipc-dashboard.org/).

## Curated data
Curated data and various mapping and translation files are placed in ./source_data.
The primary input data processed by the R script pipeline is stored in tab-delimited text files, one for each type of curated data.  These are currently:
* HIPC Dashboard - Gene Expression.tsv
* HIPC Dashboard - Cell type Frequency.tsv

The same data is also available in an Excel file, "HIPC Dashboard.xlsx", for easy review, however this file is not used by the script.

## R scripts
R scripts are in ./src
"generate_HIPC_submissions.R" is the main script.  It expects to be called from its source location.
It has a number of optional settings at the beginning.  In particular, one can choose to run the script for the gene data or the cell-type data.

## Output files
### Submission files
If it does not already exist, a directory called "../submissions" will be created relative to ./src.
The submission files are tab-delimited.
Gene submission files will be placed in a subdirectory called "../submissions/hipc_gene/".
Cell-type submissions will be placed in a subdirectory called "../submissions/hipc_ctf/".

For easy inspection, CSV-formatted versions of the same files are also generated under "../submissions/hipc_gene_csv/" and "../submissions/hipc_ctf_csv/"

### Log files
An additional directory (not checked-in to GitHub), "../logfiles" is created relative to the ./src directory.  After the script has run, this directory contains a number of files tracing the data transformations and final summary data, as well as "recreated_template" files that represent the "HIPC Dashboard.xlsx" sheets after all updates and transformations, in the same format as the original data.

### Reformatted data files
Certain files from the "logfiles" directory are checked in to GitHub in the "reformatted_data" directory.  These are the "recreated_template" files in Excel, RDS and tab-delimited formats, as well as the response components for each signature in Broad GMT format.  These files represent the data after all updates and transformations. The "recreated_template" files are in the same spreadsheet format as the original data.



