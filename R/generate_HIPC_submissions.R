# File name: generate_HIPC_submissions.R
# Author: Kenneth C. Smith

# Description:
#  Starting with curated HIPC data, prepare Dashboard submission templates.
#  Types currently explicitly supported: Gene, Cell type frequency.
#  * Perform gene symbol updates to conform to HGNC.
#    * Fix up temporary term problems pending e.g. ontology updates
#    * When exposure material is an influenza vaccine,
#        substitute the actual virus components for the listed vaccine year
#        into the target pathogen field.
#    * Substitute in corrected cell type and cell type details from external sheet.
#    * Retrieve publication title and abstract etc. based on PMID
#    * Create individual submission files by publication
#    * Create CV-per-template file for set of submissions
#
# Input files (must be in current directory):
#   * HIPC Dashboard - Gene Expression.tsv
#   * HIPC Dashboard - Cell type Frequency.tsv
#   * "vaccine_years.txt" - maps vaccine season year to vaccine viral components
#   * "manual_gene_symbol_corrections.txt" - maps invalid symbols  to known valid sybmols
#   * "cell_type_frequency-response_components_mapping.txt" - corrections to exposure_material
#        terms for cell type frequency sheet.  These corrections will in the end be added
#        directly into the HIPC Dashboard spreadsheet.
#
# Output files:
#   * Submission file for each signature
#   * Files containing the list of response_components for each signature.
#       One file has only the HGNC symbols, the other has all original symbols.
#   * CV-per-template file for all signatures of a type (gene, cell-type)
#   * Numerous log files with details of each stage of processing.
#

# BiocManager::install("xlsx")
# BiocManager::install("splitstackshape")
# BiocManager::install("HGNChelper")
# BiocManager::install("limma")
# BiocManager::install("uniqtag")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("easyPubMed")
# BiocManager::install("data.table")
# BiocManager::install("ggplot2")
# BiocManager::install("xml2")
# BiocManager::install("annotate")
# BiocManager::install("stringi")
# BiocManager::install("pracma")

library(xlsx)
library(splitstackshape)  # for cSplit()
library(uniqtag)          # for cumcount()
# Note that rbindlist() tends to return a data.table with a data.frame, which causes
# no end of problems if not cast back to just a data.frame.
library(data.table)       # for rbindlist()
library(ggplot2)
library(stringi)

source("hipc_utils.R")
source("gene_routines.R")
source("vaccines_pathogens.R")
source("pmid_to_title_easy.R")
source("write_submissions.R")
source("msigdb_submission_utils.R")
source("find_unique.R")

#####<<<< START HERE >>>>#####
##### Choose a sheet type (from "HIPC Dashboard.xlsx") #####
# Available sheet_type values are "GENE", "CELLTYPE_FREQUENCY"
sheet_type <- "CELLTYPE_FREQUENCY"

# For the moment, assume executing interactively from the ./R directory
source_data_dir <- "../source_data"
submission_dir  <- "../submissions"
logdir          <- "../logfiles"

vaccine_tsv     <- "vaccine_years.txt"
ctf_fixes_tsv   <- "cell_type_frequency-response_components_mapping.txt"
manual_gene_corrections_tsv <- "manual_gene_symbol_corrections.txt"

# symbols that have no representation in NCBI
ncbi_no_symbol <- c("BACH1-IT1")  # BACH1 Intronic Transcript 1 (uncategorized)

# NCBI not using current HGNC symbol
ncbi_fixes <- data.frame(ncbi = "TRNS1", hgnc = "MT-TS1")

##### Set runtime parameters #####
# Download publication references again using PMIDs,
#   set to FALSE to reuse existing file
#   Run this every time new publications are added to the spreadsheet,
#   for each response_component type.
RENEW_PMIDS             <- FALSE

## Please update gene files before each release
## These files will be overwritten if update is requested
# Download a new copy of NCBI gene info file
DOWNLOAD_NEW_NCBI_GENES <- FALSE
# Download a new copy of the official HGNC gene mapping
DOWNLOAD_NEW_HGNC       <- FALSE

# Generate a new copy of the mSigDB submission file
CREATE_MSIGDB           <- FALSE

# In the observation summary, do not display pathogens if the vaccine already
# uses the pathogen in its name:
# VO_0004810: 2011-2012 trivalent inactivated vaccine (A/California/7/09 (H1N1,),
#     A/Perth /16/2009 (H3N2), and B/Brisbane/60/2008).
# VO_0004899: 2012-2013 seasonal trivalent inactivated influenza vaccine
#     (A/California/7/2009 (H1N1), A/Victoria/361/2011 (H3N2),
#     and B/Wisconsin/1/2010)
# VO_0004903: Inactivated monovalent influenza A/H5N1
#     (3.75 mcg hemagglutinin [HA] A/Indonesia/05/2005) split-virus (SV) vaccine (Sanofi)
vaccine_VO_has_pathogens <- c("VO_0004810", "VO_0004899", "VO_0004903")

##############################
##### Initial data setup #####
##############################
if(!dir.exists(source_data_dir)) {
  stop(paste("source data directory not found:", source_data_dir))
}

# create log files directory
if(!dir.exists(logdir)) {
  dir.create(logdir)
}

# create submission files directory
if(!dir.exists(submission_dir)) {
  dir.create(submission_dir)
}

# Note - if running interactively, will exit script in submission_dir
setwd(submission_dir)

if (sheet_type == "GENE" && DOWNLOAD_NEW_NCBI_GENES) {
  if(!update_ncbi_homo_sapiens_gene_info(source_data_dir)) {
    print("update of NCBI gene info file failed")
  }
}

options(stringsAsFactors = FALSE)  # unfortunately doesn't help with cSplit output
summary_df <- data.frame()  # initialize summary log

##### Set up file and template name components #####
# change sheet name spaces to underscores
if (sheet_type == "GENE") {
  sheet_file     <- "HIPC Dashboard - Gene Expression.tsv"
  sheet_name     <- "Gene Expression"
  sheet_name_out <- "gene_expression"
  base_filename  <- "gene_expression"
  template_name  <- "hipc_gene"
  project        <- "Gene expression response to vaccine exposure"
} else if (sheet_type == "CELLTYPE_FREQUENCY") {
  sheet_file     <- "HIPC Dashboard - Cell type Frequency.tsv"
  sheet_name     <- "Cell type Frequency"
  sheet_name_out <- "cell_type_frequency"
  base_filename  <- "cell_type"
  template_name  <- "hipc_ctf"
  cell_mapping_sheet_name   <- "HIPC_Dashboard-Cell_type_Freque"
  project        <- "Immune cell-type frequency response to vaccine exposure"
} else {
  stop("unknown sheet type")
}
# for joint summary, list all possible values of base_filename from above
all_response_types <- c("gene_expression", "cell_type")
pmid_file <- paste(source_data_dir,
                   paste(sheet_name_out, "titles_and_dates_df.RData", sep = "_"),
                   sep = "/")

insub <- read.delim(file =  paste(source_data_dir, sheet_file, sep = "/"),
                    strip.white = TRUE,
                    stringsAsFactors = FALSE)

vaccines_by_year <- read.delim(file = paste(source_data_dir, vaccine_tsv, sep = "/"),
                               stringsAsFactors = FALSE)

# Get rid of unused columns (empty in curated data) or those not meant to appear in the Dashboard.
# target_pathogen_taxonid and tissue_type_term_id are not consistently filled in yet
#  (added for only recent data),
# No error if named column does not exist in a particular sheet.
# Note - columns "short_comment" and "process_note" are removed later after not needed anymore
remove_cols <- c("spot_check",
                 "second_spot_check",
                 "method",                  # Remove because no values in vaccine sheets yet (new column)
                 "target_pathogen_taxonid",
                 "curator_comments")
insub <- insub[!(colnames(insub) %in% remove_cols)]
length(colnames(insub))

# get rid of empty columns
insub <- insub[!(grepl("X\\.[0-9]", colnames(insub)))]

### Make any changes to headers right at the beginning,
### before separate headers and data
# Save the original response_component values
colnames(insub)[colnames(insub) == "response_component"] <- "response_component_original"

# cSplit() changes split column values to NA if first column named "X"!
# The column name of the first column is removed
# before templates are written in write_submission_template()
colnames(insub)[1] <- "donotuse"
# Add back in the columns no longer included in the curation template, 
# plus response_component to position it more towards begin
insub <- data.frame(donotuse = insub[ , 1],  # for some reason column 1 name has to be respecified
                    submission_name = "",
                    submission_date = "",
                    template_name = "",
                    response_component = "",
                    insub[ , 2:ncol(insub)],
                    stringsAsFactors = FALSE)

# Create a new column for the corrected response_component values
# and correct the headers for the response_component_original column
# Add other new columns as needed per sheet type
if (sheet_type == "GENE") {
  
  # gene-specific headers
  insub$response_component[1:6] <- insub$response_component_original[1:6]
  insub$response_component_original[1:6] <-
    c("", "label", "observed", "", "", "response component (original gene symbol)")

} else if (sheet_type == "CELLTYPE_FREQUENCY") {
  insub$cell_ontology_id  <- ""
  insub$proterm_and_extra <- ""
  insub$pro_ontology_id   <- ""
  insub$fully_qualified_response_component <- ""
  
  # cell-type specific headers
  insub$response_component[1:6] <- insub$response_component_original[1:6]
  insub$cell_ontology_id[1:6]  <-
    c("", "label", "observed", "", "", "response component cell ontology ID")
  insub$proterm_and_extra[1:6] <-
    c("", "label", "observed", "", "", "response component protein ontology and extra terms")
  insub$pro_ontology_id[1:6]   <-
    c("", "label", "observed", "", "", "response component marker protein ontology ID")
  insub$fully_qualified_response_component[1:6] <-
    c("", "label", "observed", "", "", "response component (fully qualified)")
  insub$response_component_original[1:6]        <-
    c("", "label", "observed", "", "", "response component (original curated cell type)")
}

insub$response_comp_orig_cnt  <- ""
insub$response_comp_cnt       <- ""
insub$subm_obs_id             <- ""
insub$uniq_obs_id             <- ""
insub$row_key                 <- ""

insub$submission_name[1:6]    <- c("", "label", "background", "", "", "submission name")
insub$submission_date[1:6]    <- c("", "label", "background", "", "", "submission_date")
insub$template_name[1:6]      <- c("", "label", "background", "", "", "template_name")
insub$response_comp_orig_cnt[1:6]  <- 
                                 c("", "label", "observed",   "", "", "response component (original) count")
insub$response_comp_cnt[1:6]  <- c("", "label", "observed",   "", "", "response component count")
insub$subm_obs_id[1:6]        <- c("", "label", "background", "", "", "ID of observation within a publication (PMID) and for its submission type ")
insub$uniq_obs_id[1:6]        <- c("", "label", "background", "", "", "ID of observation within its submission type")
insub$row_key[1:6]            <- c("", "label", "background", "", "", "row key")

if (sheet_type == "CELLTYPE_FREQUENCY") {
  ctf_fixes <- read.delim(file = paste(source_data_dir, ctf_fixes_tsv, sep = "/"),
                          stringsAsFactors = FALSE)
  ctf_fixes <- ctf_fixes[1:6]
  # keep only relevant content
  ctf_fixes <- subset(ctf_fixes, original_annotation != "")
  s <- strict_char_check(ctf_fixes, "\xa0")  # character 160.  Dashboard loader does not like it.
  if(!is.null(s)) {
    print(paste("for ctf_fixes, found problems in column (row numbers do not include any header)", s))
  }
}

#########################################
##### Separate header and data rows #####
#########################################
# Separate the header rows and the data rows, as the data rows will be split on several columns.
# In the original template there are 7 header lines, but here the first row is not counted,
# as it taken as column headers.
header_rows <- insub[1:6,]
df2 <- insub[7:nrow(insub),]
# Depending on how the text file is created, the sheet may have blank rows at bottom.
nrow(df2)
df2 <- df2[df2$publication_reference_id != "", ]
nrow(df2)
summary_df <- add_to_summary(summary_df, "Data rows in original sheet", nrow(df2))


# Generate observation IDs based on the PMID field value and on original row number,
# to allow e.g. response_components to be summarized by observation later.
# These are generated before any row deletions so that can refer back to original template
pmids <- df2$publication_reference_id
# ID number of observation within its submission (based on PMID)
df2$subm_obs_id  <- cumcount(pmids)
# unique observation ID (row number in original template) within this sheet
df2$uniq_obs_id  <- seq(from = 8, length.out = nrow(df2))
# save so can check later for removed observations
uniq_obs_id_orig <- df2$uniq_obs_id
df2$row_key      <- paste(pmids, df2$subm_obs_id, df2$uniq_obs_id, sep = "_")

########################################################
##### Fix temporary problems with the curated data #####
########################################################

# Remove rows marked "skip"
# FIXME - not being used anymore in HIPC Dashboard source
skip <- grepl("^skip", df2$process_note, ignore.case = TRUE)
df2 <- df2[!skip, ]

# FIXME - Most values are empty, not "Y" or "N".
#         Empty is not the same as "N".
# Currently only used for type GENE
if(sheet_type == "GENE") {
  df2$is_model[df2$is_model == ""] <- "N"
}

# "Blood" is the official UBERON name
# this should also catch the variant "whole blood cell"
w <- grep("whole blood", df2$tissue_type, ignore.case = TRUE)
df2$tissue_type[w] <- "blood"

# Write list of pathogens before any substitutions
write_unique_list(df2$target_pathogen, logdir, base_filename, "target_pathogens_before_fixes", do_split = TRUE)

# fix capitalization (should fix in source google sheet)
df2$response_behavior_type <- tolower(df2$response_behavior_type)
df2$response_behavior      <- tolower(df2$response_behavior)

# check if  comment column has more than 255 character limit
w <- sapply(df2$comments, function(x) {nchar(x) > 255}, USE.NAMES	= FALSE)
if(any(w)) {
  stop(paste("comment too long for row(s)", paste(df2$row_key[w], collapse = ", ")))
}

# check if comment column has non-UTF8 characters
w <- !stri_enc_isutf8(df2$comments)
if(any(w)) {
  stop(paste("comment has non-utf8 character in row(s)", paste(df2$row_key[w], collapse = ", ")))
}

# remove any period from end of comparison (so works in observation summaries)
df2$comparison <- trimws(df2$comparison)
if (any(grepl("\\.$", df2$comparison))) {
  df2$comparison <- sub("\\.$", "", df2$comparison)
}

# We don't expect non-ascii characters in text version of signature_source
w <- !stri_enc_isascii(df2$signature_source)
if(any(w)) {
  stop(paste("signature_source has non-ascii character in row(s)", paste(df2$row_key[w], collapse = ", ")))
}

# Summarize response behavior strings
response_behavior_strings <- sort(unique(df2$response_behavior))

write.table(response_behavior_strings,
            file = logfile_path(logdir, base_filename, "response_behavior_strings.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE)



##########################################################
##### Other global changes to template before splits #####
##########################################################

## Create a map of VO codes to text equivalents
## Both columns must have the same number of values per row
codes <- strsplit(df2$exposure_material_id,";")  # returns a list
text  <- strsplit(df2$exposure_material, ";") # returns a list
m <- mapply(function(x, y) {length(x) == length(y)}, codes, text)
if (!all(m)) {
  # FIXME - should just stop, no reason to continue
  warning(paste("exposure code vs text length mismatch for rows rowkeys:", paste(df2$row_key[!m], collapse = ", ")))
  warning(paste("removing non-matching rows"))
  df2 <- df2[m, ]
  codes <- codes[m]
  text <- text[m]
}
codes <- trimws(unlist(codes))
text <- trimws(unlist(text))

vo_code_map <- unique(data.frame(codes, text, stringsAsFactors = FALSE))
vo_code_map <- vo_code_map[order(vo_code_map$codes), ]

write.table(vo_code_map,
            file = logfile_path(logdir, base_filename, "VO_code_mapping.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE)

dupcodes <- vo_code_map[duplicated(vo_code_map$codes), ]$codes
w <- which(vo_code_map$codes %in% dupcodes)
write.table(vo_code_map[w, ],
            file = logfile_path(logdir, base_filename, "VO_code_duplicates.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE)


############################################################################
##### Data substitution and splitting (1): response_component_original #####
############################################################################

if (sheet_type == "GENE") {
  # df2$response_component_original values become factors!
  df2 <- cSplit(df2, "response_component_original", sep = ",", direction = "long")
  df2 <- cSplit(df2, "response_component_original", sep = ";", direction = "long") # ";" not used so far
  affyHits <- grep("///", df2$response_component_original)
  write.csv(df2[affyHits, c("response_component_original", "publication_reference_id", "subm_obs_id")],
            file = logfile_path(logdir, base_filename, "affyHits.csv"),
            row.names = FALSE)

  # using " /// " leaves extra slashes, regexp "[/]{3}" does not
  # FIXME - getting warnings since updated to R 4.1:
  # Warning message:
  #   In type.convert.default(unlist(x, use.names = FALSE)) :
  #   'as.is' should be specified by the caller; using TRUE
  df2 <- cSplit(df2, "response_component_original", sep = "[/]{3}", direction = "long", fixed = FALSE)
  # the curated data does in places have spaces as separators
  df2 <- cSplit(df2, "response_component_original", sep = " ", direction = "long")
} else if (sheet_type == "CELLTYPE_FREQUENCY") {
  df2 <- cSplit(df2, "response_component_original", sep = ";", direction = "long")
}

# The cSplit() calls produce a data.table of data.frame rows.  Coerce back to just data.frame.
df2 <- as.data.frame(df2)

# remove factors from df2$response component so that for genes, can alter values.
# FIXME - any way to gain better control of this?
df2$response_component_original <- as.character(df2$response_component_original)

if (sheet_type == "GENE") {
  # check for lists of response components within one PMID that have a large overlap with one another.
  # This is intended to catch the case were one set was accidentally appended to another.
  # (It happens!)
  ft <- check_response_components_overlap(df2, unique(pmids),
                                          min_intersection = 10, min_overlap_fraction = 0.75,
                                          require_different_behaviors = TRUE, max_hits = 100)
  
  file <- logfile_path(logdir, base_filename, "overlapping_signatures.txt")
  if(length(ft) > 0) {
    write.table(ft, file = file, row.names = FALSE)
  } else {
    print("no signature overlaps found...")
    if(file.exists(file)) {
      file.remove(file = file)
    }
  }
}

# create original gene symbols "signature" list after applying manual corrections
# FIXME - need better variable name than "signatures"
if (sheet_type == "GENE") {
  # Apply manual gene corrections
  rvl <- manual_gene_corrections(df2$response_component_original,
                                 paste(source_data_dir, manual_gene_corrections_tsv, sep = "/"))
  genes <- rvl$genes
  summary_df <- rbind(summary_df, rvl$summary)

  # Fix certain gene symbols containing "orf"
  genes <- fix_orf_symbols(genes)

  # Reconstruct original signatures, after splitting by various separators
  signatures <- lapply(unique(df2$uniq_obs_id), function(uniqID) {
    genes[df2$uniq_obs_id == uniqID]
  })
} else if (sheet_type == "CELLTYPE_FREQUENCY") {

  # Reconstruct original signatures, after splitting by various separators
  signatures <- lapply(unique(df2$uniq_obs_id), function(uniqID) {
    df2$response_component_original[df2$uniq_obs_id == uniqID]
  })
}

# list of original response components before main changes, with duplicates per signature removed
signatures_uniq <- lapply(signatures, unique)
# name for indexing into signatures_uniq
names(signatures_uniq) <- unique(df2$uniq_obs_id)

uids_list <- unique(df2$uniq_obs_id)
# get original count of rows for each uniq_obs_id.
uids_cnt  <- sapply(uids_list, function(x) {sum(df2$uniq_obs_id == x)})

for (i in 1:length(uids_list)) {
  # NOTE - df2$response_comp_orig_cnt is of type "character".  Be careful!
  df2$response_comp_orig_cnt[df2$uniq_obs_id == uids_list[i]] <- uids_cnt[i]
}

changed <- mapply(function(x, n) {length(x) != length(n)}, signatures, signatures_uniq)
if (length(changed) > 0) {
  print(sum(changed))
}

# if want to see that actual number of changes per signature
# tt <- mapply(function(x, n) {if(length(x) != length(n)) {
#                              print(paste(length(x), length(n), length(x) - length(n)))
#                             }}, signatures, signatures_uniq)


# find all duplicate symbols per signature
dups <- lapply(signatures, function(x) {unique(x[duplicated(x)])})
if(length(dups) > 0) {

  # assign each dup its row_key as name
  names(dups) <- sapply(unique(df2$uniq_obs_id), function(x) {
    unique(df2[df2$uniq_obs_id == x, "row_key"])
  })

  # only keep entries that have duplicates
  dup_resp_comps <- dups[sapply(dups, function(x) {length(x) > 0})]

  # Duplicate symbols can appear e.g. when probesets are translated to gene symbols
  summary_df <- add_to_summary(summary_df, "Number of signatures with duplicated items", length(dup_resp_comps))
  summary_df <- add_to_summary(summary_df, "Total number of unique duplicated items", length(unlist(dup_resp_comps)))

  # Write out list of counts of duplicated response components by signature
  write.table(paste(names(dup_resp_comps), sapply(dup_resp_comps, length), sep = ", "),
              file = logfile_path(logdir, base_filename, "response_component_duplicates_count.csv"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Write out list of duplicated response components by signature
  # Have to write in loop because varying number of elements per row
  outfile = logfile_path(logdir, base_filename, "response_component_duplicates.txt")
  if(file.exists(outfile)) file.remove(outfile)
  d <- mapply(function(x, n) {
         write.table(paste(c(x, n), collapse = "\t"),
                file = outfile,
                row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
         }, names(dup_resp_comps),  dup_resp_comps)


  # Write out sorted list of duplicated response components by signature
  outfile = logfile_path(logdir, base_filename, "response_component_duplicates_sorted.txt")
  if(file.exists(outfile)) file.remove(outfile)
  d <- mapply(function(x, n) {
         write.table(paste(c(x, sort(n)), collapse = "\t"),
                file = outfile,
                row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
         }, names(dup_resp_comps), dup_resp_comps)
}

start_cnt <- nrow(df2)

# Update original gene symbols to latest versions (NCBI/HGNC)
if (sheet_type == "GENE") {
  # starts with "genes" manually corrected list from above
  # make sure "genes" is still in synch with df2.
  length(genes) == nrow(df2)
  rvl <- update_gene_symbols(genes, logdir, base_filename,
                             source_data_dir,
                             DOWNLOAD_NEW_HGNC)
  genes_map <- rvl$genes_map
  summary_df <- rbind(summary_df, rvl$summary)
  rm(rvl)

  # w <- which(genes_map$Symbol != genes_map$alias)

  # Save the gene map for the unmatched symbols with their uniq_obs_id
  # to add them back to complete signatures later
  unmatched_symbols_map <- cbind(genes_map, uniq_obs_id = df2$uniq_obs_id)
  unmatched_symbols_map <- unmatched_symbols_map[is.na(unmatched_symbols_map$Symbol), ]

  # copy the fixed symbols back to the main data structure df2
  df2$response_component <- genes_map$Symbol

  # get rid of genes that had no valid symbol.
  df2 <- df2[!is.na(df2$response_component), ]

  # Remove genes that are OK in HGNC but fail at NCBI
  if(exists("ncbi_no_symbol")) {
    df2 <- df2[!(df2$response_component %in% ncbi_no_symbol), ]
  }

  # Fix symbols where NCBI is not using current HGNC symbol
  if(exists("ncbi_fixes")) {
    for(i in 1:nrow(ncbi_fixes)) {
      df2$response_component[df2$response_component == ncbi_fixes$hgnc[i]] <- ncbi_fixes$ncbi[i]
    }
  }

  summary_df <- add_to_summary(summary_df,
                               "Number of valid gene symbols" ,
                               length(unique(df2$response_component)))

} else if (sheet_type == "CELLTYPE_FREQUENCY") {

  # ctf_match will have NA values where no match found
  ctf_match <- match(df2$response_component_original, ctf_fixes$original_annotation)
  # Test if any rows in cell-type mapping sheet are unused.
  w <- !(1:length(ctf_fixes$original_annotation) %in% unique(ctf_match))
  if(any(w)) {
    print("not all cell-type mappings used:")
    print(sapply(ctf_fixes$original_annotation[w], function(x) {which(ctf_fixes$original_annotation == x) + 1}))
  } else {
    print("all cell-type mappings used")
  }

  if (any(is.na(ctf_match))) {
    w <- which(is.na(ctf_match))
    print(paste("row", df2$uniq_obs_id[w], ", no mapping for cell type",
                df2$response_component_original[w]))
    missing_mappings <- data.frame(row = df2$uniq_obs_id[w], label = df2$response_component_original[w])
    write.xlsx(missing_mappings,
               file = logfile_path(logdir, base_filename, "missing_mappings.xlsx"),
               sheetName = "all rows", append = FALSE, row.names = FALSE)
    write.xlsx(unique(sort(df2$response_component_original[w])),
               file = logfile_path(logdir, base_filename, "missing_mappings.xlsx"),
               sheetName = "unique", append = TRUE, row.names = FALSE)
    stop("cell type mapping not found")
  } else {
    print("all cell types matched")
  }

  failed_matches <- df2$response_component_original[is.na(ctf_match)]
  successful_matches <- df2$response_component_original[!is.na(ctf_match)]
  summary_df <- add_to_summary(summary_df, "Unique failed cell type tag matches",
                               length(unique(failed_matches)))
  summary_df <- add_to_summary(summary_df, "Unique successful cell type tag matches" ,
                               length(unique(successful_matches)))

  if (length(failed_matches) > 0) {
    write.table(unique(failed_matches),
                file = logfile_path(logdir, base_filename, "no_cell_type_match.txt"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    write.xlsx(unique(failed_matches),
               file = logfile_path(logdir, base_filename, "no_cell_type_match.xlsx"),
               row.names = FALSE, col.names = FALSE)

    write.table(sort(failed_matches_uniq),
                file = logfile_path(logdir, base_filename, "no_cell_type_match_sorted.txt"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    write.xlsx(sort(failed_matches_uniq),
               file = logfile_path(logdir, base_filename, "no_cell_type_match_sorted.xlsx"),
               row.names = FALSE, col.names = FALSE)
  }

  # combine proterm(s) and extra together.  All same length as df2.
  proterm <- ctf_fixes$PRO_term_w_intensity[ctf_match]
  extra <- ctf_fixes$extra[ctf_match]
  # vector initialized with ""
  proterm_and_extra <- vector(mode = "character", length = length(proterm))

  for (i in 1:length(proterm)) {
    if(is.na(ctf_match[i])) {
      # should not happen, already checked
      print(paste("no mapping for cell type at index", i))
      next
    }
    # insert comma between non-empty terms
    a <- c(proterm[i], extra[i])
    a <- a[a!=""]                   # returns char(0) if both empty
    a <- paste(a, collapse = ", ")  # returns "" on empty
    proterm_and_extra[i] <- ifelse(a != "", paste("&", a), "")
  }

  df2$response_component <- ctf_fixes$CL_term[ctf_match]
  df2$cell_ontology_id   <- ctf_fixes$CL_ID[ctf_match]
  df2$proterm_and_extra  <- proterm_and_extra
  df2$pro_ontology_id    <- ctf_fixes$PRO_ID[ctf_match]
  df2$fully_qualified_response_component <- ""

  # Remove any unmatched cell types (should have already been caught)
  df2 <- df2[!is.na(ctf_match), ]

  summary_df <- add_to_summary(summary_df,
                            "Unique cell types after substitutions",
                            length(unique(df2$response_component)))
  # remove trailing white space if proterm_and_extra was empty
  concatenated_cell_types <- trimws(paste(df2$response_component, df2$proterm_and_extra, sep = " "))
  df2$fully_qualified_response_component <- concatenated_cell_types

  concatenated_cell_types <- sort(unique(concatenated_cell_types))
  write.table(concatenated_cell_types,
              file = logfile_path(logdir, base_filename, "concatenated_cell_types.txt"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  summary_df <- add_to_summary(summary_df,
                            "Unique cell type + marker combinations",
                            length(concatenated_cell_types))
}

end_cnt <- nrow(df2)
summary_df <- add_to_summary(summary_df,
                             "Split 1, first reduction, rows lost to unknown response component",
                             start_cnt - end_cnt)

# Remove exact duplicates (implies duplicate copies of response component in a signature)
# This can happen when e.g. two or more original probesets get map to one symbol.
start_cnt <- end_cnt
df2 <- unique(df2)
end_cnt <- nrow(df2)

summary_df <- add_to_summary(summary_df, "Split 1, rows before removal" , start_cnt)
summary_df <- add_to_summary(summary_df, "Split 1, rows with non-unique response component removed" , start_cnt - end_cnt)
summary_df <- add_to_summary(summary_df, "Split 1, rows remaining", end_cnt)

start_cnt <- end_cnt

# get count of rows for each uniq_obs_id.
uids_list <- unique(df2$uniq_obs_id)
uids_cnt  <- sapply(uids_list, function(x) {sum(df2$uniq_obs_id == x)})

for (i in 1:length(uids_list)) {
  df2$response_comp_cnt[df2$uniq_obs_id == uids_list[i]] <- uids_cnt[i]
}

start_cnt <- end_cnt

write_unique_list(df2$response_component, logdir, base_filename, "response_component_list")

#################################################
##### Data splitting (2): exposure material #####
#################################################

# Only gene and cell type frequency have data of this type so far
if (sheet_type == "CELLTYPE_FREQUENCY" | sheet_type == "GENE") {
  write_unique_list(df2$exposure_material, logdir, base_filename, "exposure_material", do_split = TRUE)
  df2 <- cSplit(df2, "exposure_material_id", sep = ";", direction = "long")
  write_unique_list(df2$exposure_material_id, logdir, base_filename, "exposure_material_id", do_split = FALSE)
}

end_cnt <- nrow(df2)
summary_df <- add_to_summary(summary_df, "Split 2, on exposure material, rows added", end_cnt - start_cnt)
summary_df <- add_to_summary(summary_df, "Split 2, on exposure material, total rows", end_cnt)
summary_df <- add_to_summary(summary_df, "unique VO codes in exposure_material_id", length(unique(df2$exposure_material_id)))
start_cnt <- end_cnt

# The cSplit() calls produce a data.table of data.frame rows.  Coerce back to just data.frame.
df2 <- as.data.frame(df2, stringsAsFactors = FALSE)

# Substitute in actual virus components for influenza pathogens only.
for (i in 1:nrow(df2)) {
  if(df2[i, "vaccine_year"] != "" ) {
    df2[i, "target_pathogen"] <- lookup_vaccine(as.character(df2[i, "vaccine_year"]))
  }
}

# Substitute in actual pathogens for measles + DPT vaccination
# Note - the specific pathogen strains below are not specified by the VO codes used.
# These are the only case in Gene and Cell-type where there are more than a single pathogen,
# other than for influenza, so have to rewrite the pathogen entries to match split vaccine VO codes
w <- which(df2$exposure_material_id == "VO_0000654")
df2[w, "target_pathogen"] <- "Measles virus strain Edmonston-Zagreb"
w <- which(df2$exposure_material_id == "VO_0000738")
df2[w, "target_pathogen"] <- "Corynebacterium diphtheriae; Clostridium tetani; Human poliovirus 1; Human poliovirus 2; Human poliovirus 3"
is.factor(df2$exposure_material_id)


#########################################
##### Data splitting (3): Pathogens #####
#########################################

# split
df2 <- cSplit(df2, "target_pathogen", sep = ";", direction = "long")
# class(df2)  # "data.table" "data.frame"
df2 <- as.data.frame(df2)

end_cnt <- nrow(df2)
summary_df <- add_to_summary(summary_df, "Split 3, on pathogens, rows added", end_cnt - start_cnt)
summary_df <- add_to_summary(summary_df, "Split 3, on pathogens, total rows", end_cnt)
start_cnt <- end_cnt

# Write unique list of pathogens after substitutions
write_unique_list(df2$target_pathogen, logdir, base_filename, "target_pathogens_after_fixes")
summary_df <- add_to_summary(summary_df, "target_pathogens_after_fixes", length(unique(df2$target_pathogen)))

# Now substitute for viral strains missing in NCBI Taxonomy
df2$target_pathogen <- sub("Influenza A virus \\(A/Switzerland/9715293/2013\\(H3N2\\)\\)", "H3N2 subtype", df2$target_pathogen)
df2$target_pathogen <- sub("Influenza B virus \\(B/Phuket/3073/2013\\)", "Influenza B virus", df2$target_pathogen)

# lowercase first letter of each pathogen
substr(df2$target_pathogen, 1, 1) <- tolower(substr(df2$target_pathogen, 1,1))

# Write unique list of pathogens after final substitutions
write_unique_list(df2$target_pathogen, logdir, base_filename, "target_pathogens_after_final_fixes")
summary_df <- add_to_summary(summary_df, "target_pathogens_after_final_fixes", length(unique(df2$target_pathogen)))

# Write unique list of tissue types
write_unique_list(df2$tissue_type, logdir, base_filename, "tissue_type_list")
summary_df <- add_to_summary(summary_df, "tissue types", length(unique(df2$tissue_type)))

####################################################################
#### Get counts of each response component for e.g. word cloud #####
####################################################################

response_df <- unique(df2[ , c("response_component", "uniq_obs_id")])
response_df <- as.data.frame(table(response_df$response_component))
colnames(response_df) <- c("response_component", "count")
response_df <- response_df[order(response_df$count, decreasing = TRUE), ]

#Write out counts in tab-delimited format
write.table(response_df,
            file = logfile_path(logdir, base_filename, "response_component_counts.txt"),
            sep = "\t", row.names = FALSE)

# Write out counts in Excel format
write.xlsx(response_df,
           file = logfile_path(logdir, base_filename, "response_component_counts.xlsx"),
           sheetName = sheet_name, row.names = FALSE)

# Write out counts in RDS format
saveRDS(response_df,
        file = logfile_path(logdir, base_filename, "response_component_counts.RDS"))

##########################################################
####  Collect publication titles, dates and abstracts ####
##########################################################

pmids <- unique(df2$publication_reference_id)
summary_df <- add_to_summary(summary_df, "Unique PMIDs", length(pmids))
print_pub_year <- df2$publication_year[match(pmids, df2$publication_reference_id)]


# pmids <- 16571413   # for testing
# pmids <- "24336226" # for testing
# pmids <- "23594957" # for testing
if (RENEW_PMIDS || !file.exists(pmid_file)) {
  td <- mapply(pmid_to_title_easy, pmids, print_pub_year, SIMPLIFY = FALSE)
  titles_and_dates_df <- as.data.frame(rbindlist(td))
  save(titles_and_dates_df, file = pmid_file)
} else {
  load(file = pmid_file)
}

#############################################################
#### Recreate original spreadsheet with all corrections #####
#############################################################

# Collect multi-valued columns by unique observation ID (each row in original spreadsheet)
uniqIDs                    <- unique(df2$uniq_obs_id)
resp_components_cnt_df     <- data.frame()
resp_components_annotated  <- vector("list", length(uniqIDs))
resp_components_collected  <- vector("list", length(uniqIDs))
resp_components_full_sig   <- vector("list", length(uniqIDs))   # used for cell types only
recreated_template         <- vector("list", length(uniqIDs))

# Some signatures are lost entirely during cleaning
which(!(names(signatures_uniq) %in% uniqIDs))
signatures_uniq <- signatures_uniq[names(signatures_uniq) %in% uniqIDs]


# create data structures needed for mSigDB
if (sheet_type == "GENE" && CREATE_MSIGDB) {
  msigdb_empty <- msigdb_intialize()
  msigdb_list <- vector("list", length(uniqIDs))
}

for (i in 1:length(uniqIDs)) {
  df2tmp <- df2[df2$uniq_obs_id == uniqIDs[i], ]
  # Recreate a full signature in one row
  base_row <- df2tmp[1, ] # get first row for this uniqID

  response_rowname     <- paste(base_row$publication_reference_id, base_row$subm_obs_id, uniqIDs[i], sep = "_")
  response_description <- paste("PMID", base_row$publication_reference_id, "submission", base_row$subm_obs_id, "row", uniqIDs[i], sep = " ")

  # Use the full original set of response components rather than just those
  # for which a valid symbol was found.
  base_row$response_component_original <- paste(unique(df2tmp$response_component_original), collapse = "; ")

  # The "updated" (gene or cell types after fixes) response_components
  resp_components_collected[[i]] <- unique(df2tmp$response_component)  # genes or base cell types

  if (sheet_type == "GENE") {
    base_row$response_component <- paste(unique(df2tmp$response_component), collapse = "; ")
    resp_components_annotated[[i]] <- c(response_rowname, response_description, unique(df2tmp$response_component))

    # get abstract by PMID
    w <- which(titles_and_dates_df$pmid == base_row$publication_reference_id)
    if (length(w) != 1) {
      stop(paste("unexpected PMID result: uniqID = ", uniqIDs[i], ", pmid = ", base_row$publication_reference_id))
    }
    titles_and_dates_row <- titles_and_dates_df[w, ]

    if(CREATE_MSIGDB) {
      # Create mSigDB submission row
      msigdb_df <- msigdb_process_row(msigdb_empty, df2tmp)
      msigdb_list[[i]] <- msigdb_df
    }


  } else if (sheet_type == "CELLTYPE_FREQUENCY") {
    full_sig <- unique(df2tmp$fully_qualified_response_component)
    base_row$response_component    <- paste(full_sig, collapse = "; ")
    resp_components_full_sig[[i]]  <- full_sig
    resp_components_annotated[[i]] <- c(response_rowname, response_description, full_sig)
  }

  tmp <- data.frame(rowname = response_rowname,
                    pmid = base_row$publication_reference_id,
                    subm_obs_id = base_row$subm_obs_id,
                    uniq_obs_id = uniqIDs[i],
                    count = length(resp_components_collected[[i]]))
  resp_components_cnt_df <- rbind(resp_components_cnt_df, tmp)

  # Reconstitute target_pathogen and exposure_material_id
  base_row$target_pathogen   <- paste(unique(df2tmp$target_pathogen), collapse = "; ")
  base_row$exposure_material_id <- paste(unique(df2tmp$exposure_material_id), collapse = "; ")

  recreated_template[[i]] <- base_row
}
names(resp_components_collected) <- uniqIDs  # name not actually used again?
names(resp_components_full_sig)  <- uniqIDs  # name not actually used again?
names(resp_components_annotated) <- uniqIDs  # name is used for this one.

recreated_template_df <- as.data.frame(rbindlist(recreated_template))  # consolidate to a single data.frame
if(any(colnames(header_rows) != colnames(recreated_template_df))) {
  stop("mismatch between header rows and recreated_template_df rows")
}
recreated_template_df <- rbind(header_rows, recreated_template_df)

# Finish up and write out mSigDB submission
if (sheet_type == "GENE" && CREATE_MSIGDB) {
  summary_df <- write_msigdb_submission(msigdb_list, summary_df)
}

# Delete unneeded columns from the recreated template
if (sheet_type == "GENE") {
  # do nothing
} else if (sheet_type == "CELLTYPE_FREQUENCY") {
  del_cols <- c("proterm_and_extra", "fully_qualified_response_component")
  recreated_template_df <- recreated_template_df[!colnames(recreated_template_df) %in% del_cols]
}

del_cols <- c("submission_name", "submission_date", "template_name", "short_comment", "process_note")
recreated_template_df <- recreated_template_df[!colnames(recreated_template_df) %in% del_cols]

# Set that first column name back to blank
colnames(recreated_template_df)[1] <- ""
# Write out the recreated upload template in tab-delimited format
write.table(recreated_template_df,
            file = logfile_path(logdir, base_filename, "recreated_template.txt"),
            sep = "\t", row.names = FALSE)

# Write out the recreated upload template in Excel format
write.xlsx(recreated_template_df,
           file = logfile_path(logdir, base_filename, "recreated_template.xlsx"),
           sheetName = sheet_name, row.names = FALSE)

# Write out the recreated upload template in RDS format
saveRDS(recreated_template_df,
        file = logfile_path(logdir, base_filename, "recreated_template.RDS"))



########################################
##### Prepare submission templates #####
########################################

# Format URLs and handle case of multiple URLs
# Due to length restriction of 255 characters on evidence columns,
# can in practice only return at most two concatenated URLs.
# NOTE - There are special cases in this routine
# NOTE - Update code if new publication URL formats are added
df2$signature_source_url <- format_hipc_urls(df2$signature_source_url)

# do a final check for illegal characters in df2
s <- strict_char_check(df2, "\xa0")  # character 160.  Dashboard loader does not like it.
if(!is.null(s)) {
  print(paste("for df2, found problems in column (row numbers do not include any header)", s))
}
s <- strict_char_check(df2, "\n")    # embedded newlines
if(!is.null(s)) {
  print(paste("for df2, found problems in column (row numbers do not include any header)", s))
}

# Remove columns not needed for Dashboard.
del_cols <- c("exposure_material", "short_comment", "process_note")
df2 <- df2[!colnames(df2) %in% del_cols]
header_rows <- header_rows[!colnames(header_rows) %in% del_cols]


if (sheet_type == "GENE") {
  write_submission_template(df2, header_rows, template_name, titles_and_dates_df,
                            resp_components_collected, unmatched_symbols_map,
                            sheet_type, project)
} else if (sheet_type == "CELLTYPE_FREQUENCY") {
  write_submission_template(df2, header_rows, template_name, titles_and_dates_df,
                            resp_components_full_sig, unmatched_symbols_map = NULL,
                            sheet_type, project)
}

###########################################
##### Write out various summary files #####
###########################################

summary_df <- add_to_summary(summary_df, "total signatures", length(resp_components_collected))

dd0 <- sum(duplicated(resp_components_collected))
summary_df <- add_to_summary(summary_df, "duplicate signatures (response components only)", dd0)
dd0 <- length(unique(resp_components_collected))
summary_df <- add_to_summary(summary_df, "unique signatures (response components only)", dd0)

if (sheet_type == "GENE") {
  # do nothing
} else if (sheet_type == "CELLTYPE_FREQUENCY") {
  dd0 <- length(resp_components_full_sig)
  summary_df <- add_to_summary(summary_df, "signatures (including proterms)", dd0)
  dd0 <- length(unique(resp_components_full_sig))
  summary_df <- add_to_summary(summary_df, "unique signatures (including proterms)", dd0)

  outfile = logfile_path(logdir, base_filename, "full_signatures_with_proterms.txt")
  if(file.exists(outfile)) file.remove(outfile)
  d <- lapply(resp_components_full_sig,
              function(x) {write.table(paste(x, collapse = "\t"),
                          file = outfile, row.names = FALSE, col.names = FALSE,
                          quote = FALSE, append = TRUE)})

  outfile = logfile_path(logdir, base_filename, "full_signatures_with_proterms_unique.txt")
  if(file.exists(outfile)) file.remove(outfile)
  d <- lapply(unique(resp_components_full_sig),
              function(x) {write.table(paste(x, collapse = "\t"),
                          file = outfile, row.names = FALSE, col.names = FALSE,
                          quote = FALSE, append = TRUE)})
}

write.csv(resp_components_cnt_df,
          file = logfile_path(logdir, base_filename, "response_component_counts_by_row.csv"),
          row.names = FALSE)

outfile = logfile_path(logdir, base_filename, "response_components.gmt.txt")
if(file.exists(outfile)) file.remove(outfile)
d <- lapply(resp_components_annotated,
            function(x) write.table(paste(x, collapse = "\t"),
                        file = outfile, row.names = FALSE, col.names = FALSE,
                        quote = FALSE, append = TRUE))


# Convert to GMT format using same first two columns as resp_components_annotated.
# Using "c()" encloses each item or phrase in quotes.
# verify order unchanged
all(names(signatures_uniq) == names(resp_components_annotated))
signatures_uniq_gmt <- mapply(function(x, n) {
                         c(x[1], x[2], n)
                       }, resp_components_annotated, signatures_uniq)

outfile = logfile_path(logdir, base_filename, "response_components_original.gmt.txt")
if(file.exists(outfile)) file.remove(outfile)
d <- lapply(signatures_uniq_gmt,
            function(x) write.table(paste(x, collapse = "\t"),
                                    file = outfile, row.names = FALSE, col.names = FALSE,
                                    quote = FALSE, append = TRUE))


# totals by PMID
rc_cnt_by_pmid <- lapply(pmids, function(pmid) {
  vals <- resp_components_cnt_df[resp_components_cnt_df$pmid == pmid, "count"]
  return(list(pmid = pmid, sum = sum(vals)))
})
rc_cnt_by_pmid <- rbindlist(rc_cnt_by_pmid)
# Note that response components can appear in more than one signature per PMID.  These are not unique counts.
write.csv(rc_cnt_by_pmid,
          file = logfile_path(logdir, base_filename, "response_component_non-unique_count_by_PMID.csv"),
          row.names = FALSE)

# Unique response component count by PMID
# ordered by count, but could as well be ordered by PMID
rc_cnt_by_pmid_uniq <- resp_comps_by_pmid(resp_components_annotated, pmids)
rc_cnt_by_pmid_uniq <- rc_cnt_by_pmid_uniq[order(rc_cnt_by_pmid_uniq$count, decreasing = TRUE), ]

write.csv(rc_cnt_by_pmid_uniq,
          file = logfile_path(logdir, base_filename, "response_component_unique_count_by_PMID.csv"),
          row.names = FALSE)

summary_df <- add_to_summary(summary_df, "min signature size", min(resp_components_cnt_df$count))
summary_df <- add_to_summary(summary_df, "max signature size", max(resp_components_cnt_df$count))

# uniq_obs_ids that got deleted because no valid response_component
s <- subset(uniq_obs_id_orig, !(uniq_obs_id_orig %in% uniqIDs))
s <- paste(s, collapse = " ")
summary_df <- add_to_summary(summary_df, "original template rows deleted", s)

write.csv(summary_df,
          file = logfile_path(logdir, base_filename, "run_summary_counts.csv"),
          row.names = FALSE)

# if have all available response types, write a joint summary
write_final_summary(all_response_types, logdir)

##########################
##### Generate Plots #####
##########################

# Plot a histogram of unique response components per PMID
# outfile <- paste(base_filename, "UniqueCountByPMID.png", sep = "-")
# png(outfile)
# if (sheet_type == "GENE") {
#   p <- ggplot(rc_cnt_by_pmid_uniq, aes(count)) + geom_histogram(bins = 10) + scale_x_log10()
#   p +  labs(x = "genes per PMID") + annotation_logticks(sides = "b")
# } else  if (sheet_type == "CELLTYPE_FREQUENCY") {
#   p <- ggplot(rc_cnt_by_pmid_uniq, aes(count)) + geom_histogram(bins = 10)
#   p +  labs(x = "Cell types per PMID")
# }
# dev.off()

# Plot a histogram of response_components per signature
png(logfile_path(logdir, base_filename, "unique_count_by_signature.png"))
if (sheet_type == "GENE") {
  p2 <- ggplot(resp_components_cnt_df, aes(count)) + geom_histogram(bins = 10) + scale_x_log10()
  p2 +  labs(x = "genes per signature") + annotation_logticks(sides = "b") + theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
} else  if (sheet_type == "CELLTYPE_FREQUENCY") {
  p2 <- ggplot(resp_components_cnt_df, aes(count)) + geom_histogram(bins = 9) + scale_x_continuous(breaks = seq(1, 9, 1))
  p2 +  labs(x = "Cell types per signature") + theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
}
dev.off()

#################################
##### Write out sessionInfo #####
#################################

writeLines(capture.output(sessionInfo()),
           logfile_path(logdir, base_filename, "session_info.txt"))

