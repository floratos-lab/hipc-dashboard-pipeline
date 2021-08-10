# HIPC Dashboard Pipeline
## Overview
This respository, hipc-dashboard-pipeline, provides code and data to generate submission files for the HIPC Dashboard (http://hipc-dashboard.org/).  The HIPC Dashboard itself resides in a separate GitHub repository at https://github.com/floratos-lab/hipc-signature.

The HIPC Dashboard provides a web interface to the immune signatures curated as part of the HIPC Signatures II project (NIAID).  This initial version of the Dashboard focuses on vaccine reponse signatures, however, it is designed to be extendable to additional signature types.  Work on signatures of responses to infection is in progress.  

The focus of curation are the "response components", the biological features whose change is being measured.  Extensive supporting metadata is captured to characterize the experiment giving rise to each signtaure.  The intial response components curated and released on the Dashboard have been gene expression and cell-type frequency.  Work is also in progress on metabolites, pathways etc.

## Design
The Dashboard design is flexible and allows submissions with an arbitrary number and type of data columns.  The Dashboard is built on two classes of data.  The first, termed "subjects", comprises terms drawn from controlled vocabularies which are represented directly in the Dashboard database with all their underlying data. The second class, termed "evidence", is open and can be used to add additionnal annotation columns, such as free text or files, as needed for a particular submission.   The addition of new ontology-based data-types ("subjects") requires additions to the data model.  The templates for the first two supported submission types, gene expression and cell-type frequency, are almost identical. The HIPC Dashboard design is based on an earlier project, the Cancer Target Detection and Discovery Network (CTD2), initiated by the Office of Cancer Genomics of the National Cancer Institute (NCI).  A more detailed discussion of the Dashboard architecture is available in the publication Askoy et al. (2017), https://pubmed.ncbi.nlm.nih.gov/29220450/.

## Curation
In general, curators initially enter data into type-specific templates exactly as it appears in a publication.  This data is then standardized as needed using suitable ontologies.  The original annotations are preserved for qualtity control and provenance.  The curation sheets contain several rows of headers used to guide the Dashboard load process but which are not of interest to the wider community.

## Data Standardization
Standardization varies according to data type.  We use existing community standards wherever possible:
* Gene symbols - Curated gene symbols are updated to current HGNC/NCBI symbols based on (1) NCBI synonyms and (2) a manually created mapping table.  The later deals with specific problematic symbols found in the curated data, where examination of the original data is able to support a definite mapping.
* Cell types and markers - We create a mapping table to standardize the orignal cell type descriptions using terms from the Cell Ontology for cell types and Protein Ontology for additional type-defining markers.
* Vaccines - For influenza vaccines, the year is used to expand the vaccine into its three or four viral components.

## Curated data repository
Curated data and various mapping and translation files are placed in ./source_data.
The primary input data processed by the R script pipeline is stored in tab-delimited text files, one for each type of curated data.  These are currently:
* HIPC Dashboard - Gene Expression.tsv
* HIPC Dashboard - Cell type Frequency.tsv

The same data is also available in an Excel file, "HIPC Dashboard.xlsx", for easy review, however this file is not used by the script.

## R scripts
R scripts are in ./src.
"generate_HIPC_submissions.R" is the main script.  It expects to be called from its source location.
It has a number of optional settings at the beginning.  In particular, one can choose to run the script for the gene data or the cell-type data.

## Data Releases
The Dashboard database is reloaded in its entirety each time a new release is created.  When a new version of the data is ready, the code, input data and resulting output files are all committed together as a release.  

## Output files
### Submission files (Dashboard load files)
If it does not already exist, a directory called "../submissions" will be created relative to ./src.
The submission files are tab-delimited.
Gene submission files are placed in a subdirectory called "../submissions/hipc_gene/".
Cell-type submissions are placed in a subdirectory called "../submissions/hipc_ctf/".

For easy inspection, CSV-formatted versions of the same files are also generated under "../submissions/hipc_gene_csv/" and "../submissions/hipc_ctf_csv/"

### Log files
An additional directory (not checked-in to GitHub), "../logfiles" is created relative to the ./src directory.  After the script has run, this directory contains a number of files tracing the data transformations and final summary data, as well as "recreated_template" files that represent the data after all updates and transformations, in the same format as the original data.

### Reformatted data files
The "recreated_template" files are in the same spreadsheet format as the original curated data.  These files represent the original data after all updates and transformations, and are provided in tab-delimited, RDS and Excel formats.  Files containing the "response components" for each signature are also provided in the tab-delimited Broad GMT format.

## Database Schema
![Dashboard DB schema](hipc_dashboard_data_model_simple.png)

## References
Aksoy BA, Dancík V, Smith K, Mazerik JN, Ji Z, Gross B, Nikolova O, Jaber N, Califano A, Schreiber SL, Gerhard DS, Hermida LC, Jagu S, Sander C, Floratos A, Clemons PA. CTD2 Dashboard: a searchable web interface to connect validated results from the Cancer Target Discovery and Development Network. Database (Oxford). 2017 Jan 1;2017:bax054. doi: 10.1093/database/bax054. PMID: 29220450; PMCID: PMC5569694.
