
The processing pipeline consists of four scripts `main_inf_ctf.R`, `main_inf_gene.R`, `main_vac_ctf.R`, and `main_vac_ctf.R`.
The input data for the scripts are in the directory `data/source_curations` and `data/references`; 
the output goes to `data/release_files`, which will be used by the dashboard application data loading process.

The four scripts are for two exposure types: infection and vaccine, 
combined with two response types: cell type frequency and gene expression.
There are four corresponding subdirectories under `data/release_files` for the results.

The files under `data/references` should be updated separately before running the main scripts.
