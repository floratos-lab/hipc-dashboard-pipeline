# depends on "source("hipc_utils.R")"  having been called in main routine
library(HGNChelper)
library(org.Hs.eg.db)
library(limma)
library(annotate)
library(RCurl)

ncbi_gene_file_ftp <- "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
hgnc_file <- "hgnc_complete_set.RData"

update_ncbi_homo_sapiens_gene_info <- function(source_data_dir) {
  if(url.exists(ncbi_gene_file_ftp)) {
    download.file("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
                  dest = paste(source_data_dir, basename(ncbi_gene_file_ftp), sep = "/"), mode = "wb")
  } else {
    return(FALSE)
  }
  return(TRUE)
}

manual_gene_corrections <- function(genes, corrections_file) {
  summary_df <- data.frame()  # initialize summary log
  # First make some manual changes
  manual_genes <- read.table(file = corrections_file, header = TRUE,
                             strip.white = TRUE,
                             stringsAsFactors = FALSE, sep = "\t")

  # How many of the correct manual symbols are already present?
  m <- match(genes, manual_genes$symbol)
  w <- which(!is.na(m))  # w shows which rows of df2 to fix
  lenPre <- length(genes[w])

  # Now change the incorrect versions
  m <- match(genes, manual_genes$alias)
  w <- which(!is.na(m))
  genes[w] <- manual_genes$symbol[m[w]]

  # Now how many of the correct manual symbols are present?
  w <- which(!is.na(match(genes, manual_genes$symbol)))  # w shows which rows of df2 to fix
  lenPost <- length(genes[w])
  summary_df <- add_to_summary(summary_df,
                            "Manual gene symbol substitutions: symbols in list", nrow(manual_genes))
  summary_df <- add_to_summary(summary_df,
                            "Manual gene symbol substitutions", lenPost - lenPre)

  return(list(genes = genes, summary = summary_df))
}

fix_orf_symbols <- function(genes) {
  # Change "ORF" to "orf"  except for genes that have ORF in them but are not orf genes,
  # e.g. MORF4L1, "Mortality Factor 4 Like 1"
  w <- grep("ORF", genes)
  skip <- grep("MORF4L", genes)
  w <- w[!(w %in% skip)]
  genes[w] <- sub("ORF", "orf", genes[w])
  return(genes)
}

# Check gene symbol list against real NCBI symbols
# Returns index of non-matching symbols
check_against_real_ncbi <- function(gene_list, source_data_dir) {
  summary_df <- data.frame()  # initialize summary log

  Homo_sapiens.gene_info.file <- paste(source_data_dir, basename(ncbi_gene_file_ftp), sep = "/")

  ncbi_genes <- read.table(file = Homo_sapiens.gene_info.file, header = TRUE,
                           strip.white = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "", comment.char = "")

  return(which(!(gene_list %in% ncbi_genes$Symbol)))
}


update_symbols_ncbi <- function(genes, source_data_dir, logdir, base_filename) {
  summary_df <- data.frame()  # initialize summary log


  Homo_sapiens.gene_info.file <- paste(source_data_dir, basename(ncbi_gene_file_ftp), sep = "/")

  # treat "genes" as aliases, find official symbols
  genes_map <- alias2SymbolUsingNCBI(genes,
                                     gene.info.file = Homo_sapiens.gene_info.file,
                                     required.columns =
                                       c("GeneID", "Symbol", "Synonyms",
                                         "Symbol_from_nomenclature_authority"))

  genes_map <- as.data.frame(cbind(genes, genes_map), stringsAsFactors = FALSE)
  names(genes_map)[1] <- "alias"

  # Aliases for which a valid symbol was found
  # Note that genes_map$Symbol has many NAs but those rows do not
  # emerge in the subset() output.  They would appear if
  # the selection was done directly.
  changedSymbols <- unique(subset(genes_map, alias != Symbol))
  changedSymbols <- changedSymbols[order(changedSymbols$alias), ]
  summary_df <- add_to_summary(summary_df,
                               "alias2SymbolUsingNCBI(): symbols changed" , nrow(changedSymbols))

  outFile <- paste0(logdir, "/", base_filename, "-bioc_substitute_genes.txt")
  write.table(changedSymbols, outFile, sep = "\t", row.names = FALSE, quote = FALSE)
  outFile <- paste0(logdir, "/", base_filename, "-bioc_substitute_genes.xlsx")
  if (nrow(changedSymbols) > 0) {
    write.xlsx(changedSymbols, outFile, row.names = FALSE)
  }

  # if no symbol was found for an alias, it is NA
  failed_symbols <- unique(subset(genes_map, is.na(Symbol)))
  failed_symbols <- sort(failed_symbols$alias)
  summary_df <- add_to_summary(summary_df,
                               "alias2SymbolUsingNCBI(): failed symbols remaining" , length(failed_symbols))

  # Are all non-matching symbols NA?
  w <- check_against_real_ncbi(genes_map$Symbol, source_data_dir)

  summary_df <- add_to_summary(summary_df,
                               "update_symbols_using_ncbi: all non-matching symbols are NA",
                               all(is.na(genes_map$Symbol[w])))

  # there are NCBI symbols that do not have HGCN symbols, which are reported as "-"
  # FIXME - should handle case where all match
  w <- which(genes_map$Symbol != genes_map$Symbol_from_nomenclature_authority)
  ncbi_vs_hgnc <- data.frame(ncbi = genes_map$Symbol[w],
                             hgnc = genes_map$Symbol_from_nomenclature_authority[w],
                             stringsAsFactors = FALSE)

  ncbi_vs_hgnc <- unique(ncbi_vs_hgnc)
  if(!all(ncbi_vs_hgnc$hgnc == "-")) {
    print(paste("update_symbols_using_ncbi: ncbi_vs_hgnc not \"-\""))
  }
  summary_df <- add_to_summary(summary_df,
                               "update_symbols_using_ncbi: number of ncbi vs hgnc differences",
                               nrow(ncbi_vs_hgnc))
  summary_df <- add_to_summary(summary_df,
                               "update_symbols_using_ncbi: all ncbi vs hgnc differences are \"-\"",
                              all(ncbi_vs_hgnc$hgnc == "-"))

  write.table(failed_symbols,
              file = logfile_path(logdir, base_filename, "bioc_no_substitute_genes.txt"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  return(list(genes_map = genes_map, failed_symbols = failed_symbols, summary = summary_df))
}

update_symbols_using_hgnc_helper <- function(genes_map, source_data_dir, logdir, base_filename, download_new_hgnc) {
  summary_df <- data.frame()  # initialize summary log
  hgnc_file_path <- paste(source_data_dir, hgnc_file, sep = "/")
  if (download_new_hgnc) {
    # Data source is "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"
    new.hgnc.table <- getCurrentHumanMap()
    save(new.hgnc.table, file = hgnc_file_path, compress="bzip2")
  }
  if(file.exists(hgnc_file_path)) {
    load(hgnc_file_path)  # assumes it was previously downloaded
    genes.retried <- checkGeneSymbols(genes_map$alias, map=new.hgnc.table, species="human")
  } else {
    print(paste("hgnc file not found:", hgnc_file_path))
    print(paste("hgnc file not found, using HGNCHelper package version"))

    # use the map that is included with the HGNChelper package.
    genes.retried <- checkGeneSymbols(genes_map$alias, species="human")
  }

  newHit <- is.na(genes_map$Symbol) & !is.na(genes.retried$Suggested.Symbol)
  # not going to try to fix these (just one present)
  newHit <- newHit & !grepl(" /// ", genes.retried$Suggested.Symbol)

  additional <- unique(genes.retried[newHit, ])
  additional <- additional[order(additional$x), ]
  summary_df  <- add_to_summary(summary_df, "HGNChelper: symbols changed" , nrow(additional))

  write.table(additional,
              file = logfile_path(logdir, base_filename, "HGNC_additional_substitute_genes.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.xlsx(additional,
             file = logfile_path(logdir, base_filename, "HGNC_additional_substitute_genes.xlsx"),
             row.names = FALSE)

  # Add the HGNC symbols to the master list
  genes_map$Symbol[newHit] <- genes.retried$Suggested.Symbol[newHit]

  # find HGNC symbols not in NCBI
  w <- check_against_real_ncbi(additional$Suggested.Symbol, source_data_dir)
  summary_df  <- add_to_summary(summary_df,
                                "HGNChelper: HGNC symbols not in NCBI" ,
                                paste(additional$Suggested.Symbol[w], collapse = ", "))

  # gene symbols still not found after HGNC helper
  failed_symbols <- subset(genes_map, is.na(Symbol))
  failed_symbols <- sort(unique(failed_symbols$alias))
  summary_df <- add_to_summary(summary_df,
                               "HGNChelper: failed symbols remaining" , length(failed_symbols))

  # Here we just look at disagreements between the bioconductor symbols and the HGNC symbols in general.
  # This does not feed into the actual final set of gene symbols we will use.
  alternate <- genes_map$Symbol != genes.retried$Suggested.Symbol
  alternate_symbols <- data.frame(symbol = genes_map$Symbol[alternate],
                                  hgnc = genes.retried$Suggested.Symbol[alternate])
  alternate_symbols <- unique(subset(alternate_symbols, !is.na(symbol)))
  summary_df <- add_to_summary(summary_df,
                               "HGNC_vs_bioc: count of differences" , nrow(alternate_symbols))

  write.table(alternate_symbols,
              file = logfile_path(logdir, base_filename, "HGNC_vs_bioc_alternate_genes.txt"),
              sep = "\t", row.names = FALSE)
  write.xlsx(alternate_symbols,
             file = logfile_path(logdir, base_filename, "HGNC_vs_bioc_alternate_genes.xlsx"),
             row.names = FALSE, col.names = FALSE)

  return(list(genes_map = genes_map, failed_symbols = failed_symbols, summary = summary_df))
}

# test_symbols <- c("AARS", "ABP1", "ACCN2", "ACPL2", "ACPP", "ACRC", "ADC", "ADCK3", "ADFP")

#### Two routines to update gene symbols by converting them to EntrezIDs and
#### then back to gene symbols using org.Hs.eg.db data

# Note - calling unlist() on a list with NULL values will drop those items,
# so list size will change.  So don't do it.  Change them to NAs.
# NOTE - org.Hs.eg.db data is not updated regularly.  It is very out-of-date at present.
# Returns same number of rows as in input

symbols_to_entrezid_org.Hs <- function(qsymbols) {
  summary_df <- data.frame()  # initialize summary log

  #### Convert to EntrezIDs
  x <- org.Hs.egSYMBOL2EG
  mapped_symbols <- mappedkeys(x)
  # Convert to a list
  # names are symbols, values are entrezids
  symb_entrezid_list <- as.list(x[mapped_symbols])  # names are symbols, values are entrezids

  x <- org.Hs.egALIAS2EG
  mapped_aliases <- mappedkeys(x)
  # Convert to a list
  # names are symbols, values are entrezids
  alias_entrezid_list <- as.list(x[mapped_aliases])

  # Change known symbols to EntrezID
  # names are gene symbols, values are entrez ids
  # unmatched list item names are <NA>, list item values are NULL
  symb_entrez_ids <- symb_entrezid_list[qsymbols]
  symb_entrez_ids <- sapply(symb_entrez_ids, "[", 1)
  # Change NULL values to NA
  w <- is.na(names(symb_entrez_ids))
  symb_entrez_ids[w] <- NA
  # now safe to call unlist()
  symb_entrez_ids <- unlist(symb_entrez_ids)
  if (length(symb_entrez_ids) != length(qsymbols)) {
    stop("symbols_to_entrezid_using_org.Hs.eg.db_preserve: symbols mapping length failure")
  }

  # Attempt to match aliases to EntrezIDs
  # names are gene symbols, values are entrez ids
  # unmatched list item names are <NA>, list item values are NULL
  alias_entrez_ids <- alias_entrezid_list[qsymbols]
  alias_entrez_ids <- sapply(alias_entrez_ids, "[", 1)
  # Change NULL values to NA as follows
  w <- is.na(names(alias_entrez_ids))
  alias_entrez_ids[w] <- NA
  # now safe to call unlist()
  alias_entrez_ids <- unlist(alias_entrez_ids)
  if (length(alias_entrez_ids) != length(qsymbols)) {
    stop("symbols_to_entrezid_using_org.Hs.eg.db_preserve: aliases mapping length failure")
  }

  #### Now weave the two results back together
  w <- is.na(symb_entrez_ids)
  symb_entrez_ids[w] <- alias_entrez_ids[w]

  na_values <- sum(is.na(symb_entrez_ids))
  summary_df <- add_to_summary(summary_df,
                             "  update_symbols_using_org.Hs.eg.db: symbols with no matching entrez ids",
                             na_values)
  return(symb_entrez_ids)
}

# NOTE - org.Hs.eg.db data is not updated regularly.  It is very out-of-date at present.
# Returns same number of rows as in input
entrezids_to_symbol_org.Hs <- function(qentrezid) {

  #### Convert entrez id to official gene symbol
  x <- org.Hs.egSYMBOL
  symbol_keys <- mappedkeys(x)
  # Convert to a list
  # Names are EntrezIDs, values are gene symbols
  symbol_map <- as.list(x[symbol_keys])
  symbols <- symbol_map[qentrezid]
  sl <- sapply(symbols, length)
  if(any(!is.na(sl) && sl != 1)) {
    stop("entrezids_to_symbol_using_org.Hs.eg.db_preserve - unexpected multiple mapping")
  }
  symbols[is.na(names(symbols))] <- NA
  symbols <- unlist(symbols)
  if (length(symbols) != length(qentrezid)) {
    stop("entrezids_to_symbol_using_org.Hs.eg.db_preserve: EntrezID mapping length failure")
  }
  return(symbols)
}

genbank_acc_to_ncbi <- function(genes_map, failed_symbols, source_data_dir) {
  ##### Handle Genbank Accession IDs like "AB004857"
  summary_df <- data.frame()  # initialize summary log

  x <- org.Hs.egSYMBOL
  # Get the gene symbols that are mapped to entrez gene identifiers
  symbol_key <- mappedkeys(x)
  # Convert to a list
  symbol_map <- as.list(x[symbol_key])

  # Change matching Genbank accession IDs to EntrezID and then look up NCBI/HGNC symbol
  # FIXME how come not doing mappedkeys(x) step here?
  accToEntrez <- as.list(org.Hs.egACCNUM2EG)
  m       <- match(failed_symbols, names(accToEntrez))  # note - there are refseq ids in here, e.g. "NM_001085"
  entrez  <- unlist(accToEntrez[m])  # the accession is the name of the entrez id value.
  symbol2 <- symbol_map[entrez]  # via entrez IDs
  acc_df   <- data.frame(acc = names(entrez), entrezid = as.integer(entrez), ncbi = as.character(symbol2))
  summary_df <- add_to_summary(summary_df,
                               "genbank: unique accessions fixed" , length(symbol2))

  # Which PMIDs do the genbank accession IDs come from
  # m <- match(names(entrez), df2$response_component)
  # print(table(df2$publication_reference[m]))
  # 16571413 26726811
  #  64        1
  # which genbank accession is the single hit?
  # w <- which(df2$publication_reference[m] == "26726811")
  # (df2$response_component[m])[w]  # "BC017937"

  # Check that genbank accession mapping got same symbols we would have gotten in the regular pipeline above..
  # new_symbols <- alias2SymbolTable(acc_df$hgnc, species = "Hs")
  Homo_sapiens.gene_info.file <- paste(source_data_dir, basename(ncbi_gene_file_ftp), sep = "/")

  new_symbols <- alias2SymbolUsingNCBI(acc_df$ncbi,
                                      gene.info.file = Homo_sapiens.gene_info.file,
                                      required.columns = c("Symbol", "Symbol_from_nomenclature_authority"))

  w <- which(new_symbols$Symbol != acc_df$ncbi)
  summary_df <- add_to_summary(summary_df,
                              "genbank: count of corrected accessions that did not get expected NCBI symbols",
                              sum(new_symbols$Symbol != acc_df$ncbi))
  if(length(w) > 0 && w > 0) {
    summary_df <- add_to_summary(summary_df,
                                "genbank: incorrect accession symbols",
                                acc_df$ncbi[w])
    summary_df <- add_to_summary(summary_df,
                                "genbank: corrected accession symbols",
                                new_symbols$Symbol[w])
    # fix any discrepancies
    acc_df$ncbi[w] <- new_symbols$Symbol[w]
  }
  # Now change the incorrect versions
  m <- match(genes_map$alias, acc_df$acc)
  w <- which(!is.na(m))  # w shows which rows of genes_map to fix
  # if want to look at the mapping
  # dd <- data.frame(alias = genes_map$alias[w], symbol = acc_df$ncbi[m[w]])
  genes_map$Symbol[w] <- acc_df$ncbi[m[w]]

  summary_df <- add_to_summary(summary_df,
                               "genbank: total accessions converted to NCBI" , length(w))

  # Example - check that the substitution worked for "Y13115"
  # print(subset(genes_map, genes_map$alias == "Y13115")) # should show the alias was replaced

  return(list(genes_map = genes_map, summary = summary_df))
}

check_against_ncbi_synonyms <- function(failed_symbols, source_data_dir, logdir, base_filename) {

  ## Check against the actual NCBI background data.
  # Here we will use grep rather than exact matches, to see if any
  # genes warrant manual investigation.
  # Note we don't use this source as primary because the bioconductor method provides
  # a defined way to chose among multiple hits.
  # "Homo_sapiens.gene_info" has "'" and "#" characters embedded in strings, e.g. "5'";
  # have to turn off "quote" and "comment.char" options to correctly read in this file.
  summary_df <- data.frame()  # initialize summary log
  Homo_sapiens.gene_info.file <- paste(source_data_dir, basename(ncbi_gene_file_ftp), sep = "/")

  ncbi_genes <- read.table(file = Homo_sapiens.gene_info.file, header = TRUE,
                           strip.white = TRUE, stringsAsFactors = FALSE, sep = "\t",
                           quote = "", comment.char = "")

  m <- match(failed_symbols, ncbi_genes$Symbol)
  sum(is.na(m))   # 929
  sum(!is.na(m))  # 0
  # Using grepl because synonyms can have many entries per line, just need to find candidates
  # Candidates can then be identified manually
  m <- sapply(failed_symbols, function(x) {any(grepl(x,ncbi_genes$Synonyms))})

  # m_lines will contain a list of vectors of matching row numbers
  # sapply keeps list item names, lapply apparently does not
  m_lines <- sapply(failed_symbols[m], function(x) {grep(x, ncbi_genes$Synonyms)})

  # process each item in a list of vectors of varying length
  print_df <- data.frame()
  cnt <- 0
  for (k in 1:length(m_lines)) {
    for (l in 1:length(m_lines[[k]])) {
      cnt <- cnt + 1
      row_n <- m_lines[[k]][l]
      # print(paste0(names(m_lines[k]), ": ",  ncbi_genes$Synonyms[row_n]))
      print_df <- rbind(print_df, data.frame(gene = names(m_lines[k]), synonyms = ncbi_genes$Synonyms[row_n]))
    }
  }

  # Write out the list of potential synonyms
  write.table(print_df,
              file = logfile_path(logdir, base_filename, "genbank_synonyms.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.xlsx(print_df,
             file = logfile_path(logdir, base_filename, "genbank_synonyms.xlsx"),
             row.names = FALSE, col.names = TRUE)
  summary_df <- add_to_summary(summary_df,
                               "ncbi: possible synonyms to be checked manually" , cnt)
  return(summary_df)
}

log_no_valid_symbol_vs_pmid <- function(noValidSymbols_df, base_filename) {
  # contains the following columns: c("response_component", "publication_reference", "subObsID", "uniqObsID")
  write.table(noValidSymbols_df,
              file = logfile_path(logdir, base_filename, "no_substitute_genes_with_PMIDs.txt"),
              row.names = FALSE, col.names = TRUE,  quote = FALSE)
  write.xlsx(noValidSymbols_df,
             file = logfile_path(logdir, base_filename, "no_substitute_genes_with_PMIDs.xlsx"),
             row.names = FALSE, col.names = TRUE)

  pmids_with_invalid_gene_symbols <- unique(noValidSymbols_df$publication_reference)
  write.table(pmids_with_invalid_gene_symbols,
              file = logfile_path(logdir, base_filename, "PMIDs_with_invalid_gene_symbols.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# Run all symbol updates
update_gene_symbols <- function(genes, logdir, base_filename, source_data_dir, download_new_hgnc) {
  summary_df <- data.frame()  # initialize summary log

  rvl <- update_symbols_ncbi(genes, source_data_dir, logdir, base_filename)
  genes_map      <- rvl$genes_map
  failed_symbols <- rvl$failed_symbols
  summary_df     <- rbind(summary_df, rvl$summary)
  rm(rvl)

  #  Warning messages for   update_symbols_using_hgnc_helper():
  #  1: In checkGeneSymbols(genes_map$alias, map = new.hgnc.table, species = "human") :
  #  Human gene symbols should be all upper-case except for the 'orf' in open reading frames. The case of some letters was corrected.
  #  2: In checkGeneSymbols(genes_map$alias, map = new.hgnc.table, species = "human") :
  #  x contains non-approved gene symbols

  # See how many NAs can be fixed
  # For some reason, HGNChelper can map aliases to HGNC symbols that ARE also NCBI symbols
  # but that alias2SymbolUsingNCBI() does not map.
  rvl            <- update_symbols_using_hgnc_helper(genes_map, source_data_dir, logdir, base_filename,
                                                     download_new_hgnc)
  genes_map      <- rvl$genes_map
  failed_symbols <- rvl$failed_symbols
  summary_df     <- rbind(summary_df, rvl$summary)
  rm(rvl)

  #### Genbank accession substitution
  rvl        <- genbank_acc_to_ncbi(genes_map, failed_symbols, source_data_dir)
  genes_map  <- rvl$genes_map
  summary_df <- rbind(summary_df, rvl$summary)

  # Final set of gene symbols not found
  failed_symbols <- subset(genes_map, is.na(Symbol))
  failed_symbols <- sort(unique(failed_symbols$alias))
  summary_df     <- add_to_summary(summary_df,
                               "Number of invalid gene symbols after all replacements", length(failed_symbols))

  # Does a check against NCBI synonyms but does not alter any values
  # Writes results to file
  summary_rv <- check_against_ncbi_synonyms(failed_symbols, source_data_dir, logdir, base_filename)
  summary_df <- rbind(summary_df, summary_rv)

  # Save PMID info for unmapped symbols
  no_valid_symbols_df <- df2[is.na(genes_map$Symbol),
                             c("response_component", "publication_reference", "subm_obs_id", "uniq_obs_id")]
  log_no_valid_symbol_vs_pmid(no_valid_symbols_df, base_filename)
  return(list(genes_map = genes_map, summary = summary_df))
}

