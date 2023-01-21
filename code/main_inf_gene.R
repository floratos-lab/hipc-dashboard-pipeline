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

library(splitstackshape)  # for cSplit()
library(uniqtag)          # for cumcount()
# Note that rbindlist() tends to return a data.table with a data.frame, which causes
# no end of problems if not cast back to just a data.frame.
library(data.table)       # for rbindlist()
library(stringr)
library(R.utils)

source("hipc_utils.R")
source("gene_routines.R")
source("pmid_to_title_easy.R")
source("write_submissions.R")

# Assume executing from the ./code directory
source_curations       <- "../data/source_curations"
reference_files        <- "../data/reference_files"

standardized_curations <- "../data/standardized_curations"
log_files              <- "../logfiles"
csv_submission_files   <- "../logfiles"

manual_gene_corrections_txt <- "manual_gene_symbol_corrections.txt"

##### Set runtime parameters #####
# Download publication references again using PMIDs,
#   set to FALSE to reuse existing file
#   Run this every time new publications are added to the spreadsheet,
#   for each response_component type.
RENEW_PMIDS             <- FALSE  

# Download a new copy of the official HGNC gene mapping
DOWNLOAD_NEW_HGNC       <- FALSE

exclude_pmid <- "33361761"  # PMID(s) that are not suitable for processing

##### Set up file and template name components #####
first_data_row <- 9  # not counting the first row, which becomes a column header.

sheet_file     <- "COVID-19 curation template - example curation.tsv"
sheet_file2    <- "Odak_2020_pmid_32650275 - covid19.tsv"
sheet_file3    <- "hipc_infection_covid_v2 - multiple_types.tsv"

base_filename  <- "inf_gene_expression"

# used to filter sheet rows
response_behavior_type_var <- "gene expression"

pmid_file <- paste(reference_files,
                   paste(base_filename, "titles_and_dates_df.RData", sep = "-"),
                   sep = "/")

insub <- read.delim(file =  paste(source_curations, sheet_file, sep = "/"),
                    strip.white = TRUE,
                    stringsAsFactors = FALSE)

# Get rid of unused columns (empty in curated data) or those not meant to appear in the Dashboard.
# No error if named column does not exist in a particular sheet.
del_cols_common <- c("submission_name", "template_name", # should not appear anymore
              "spot_check", "Spot.check",
              "second_spot_check",
              "Curator",
              "curator_comments")
insub <- insub[!(colnames(insub) %in% del_cols_common)]

# read two other source curation files
insub2 <- read.delim(file =  paste(source_curations, sheet_file2, sep = "/"),
                       strip.white = TRUE,
                       stringsAsFactors = FALSE)
insub2 <- insub2[!(colnames(insub2) %in% del_cols_common)]

insub3 <- read.delim(file =  paste(source_curations, sheet_file3, sep = "/"),
                         strip.white = TRUE,
                         stringsAsFactors = FALSE,
                         quote = "")
insub3 <- insub3[!(colnames(insub3) %in% del_cols_common)]

del_cols <- c("route",
                "addntl_time_point_units",
                "signature_source_url")  # discontinued in latest template

insub2 <- insub2[!(colnames(insub2) %in% del_cols)]
insub3 <- insub3[!(colnames(insub3) %in% del_cols)]

if( !all(colnames(insub2) == colnames(insub3)) || !all(colnames(insub) == colnames(insub2)) ) {
  colnames(insub)
  colnames(insub2)
  colnames(insub3)
  stop("column names do not match in the three source duration files")
}

# stitch together the data sections
insub <- rbind(insub, insub2[first_data_row:nrow(insub2),])
insub <- rbind(insub, insub3[first_data_row:nrow(insub3),])

### Make any changes to headers right at the beginning,
### before separating headers and data
# Save the original response_component values
colnames(insub)[colnames(insub) == "response_component"] <- "response_component_original"

# cSplit() changes split column values to NA if first column named "X"!
# The column name of the first column is removed
# before templates are written in write_submission_template()

# Add back in the columns no longer included in the curation template,
# plus response_component to position it more towards begin
insub <- data.frame(donotuse = insub[ , 1],
                    submission_name = "",
                    submission_date = "",
                    template_name = "",
                    response_component = "",
                    insub[ , 2:ncol(insub)],
                    response_comp_orig_cnt = "",
                    response_comp_cnt = "",
                    sig_subm_id = "",
                    sig_row_id = "",
                    row_key = "",
                    stringsAsFactors = FALSE)

# gene-specific headers
COLUMN_NAMES = c("donotuse", "submission_name", "submission_date", "template_name", "response_component",
  "curation_date", "cohort", "age_min", "age_max", "age_units", "number_subjects", "tissue_type", 
  "tissue_type_term_id", "method", "response_component_original", "is_model", "response_behavior_type", 
  "response_behavior", "comparison", "baseline_time_event", "time_point", "time_point_units", 
  "exposure_material_id", "exposure_process", "disease_name", "disease_stage", "publication_reference_id", 
  "publication_date", "publication_reference_url", "signature_source", "comments", 
  "response_comp_orig_cnt", "response_comp_cnt", "sig_subm_id", "sig_row_id", "row_key")
header_rows = insub[colnames(insub) %in% COLUMN_NAMES][1:6,]

# set header values in 10 columns
header_rows$response_component      <- c("gene", "", "gene_biomarker", "", "", "response component")
header_rows$response_component_original <- c("", "label", "observed", "", "", "response component (original gene symbol)")
header_rows$submission_name         <- c("", "label", "background", "", "", "submission name")
header_rows$submission_date         <- c("", "label", "background", "", "", "submission_date")
header_rows$template_name           <- c("", "label", "background", "", "", "template_name")
header_rows$response_comp_orig_cnt  <- c("", "label", "observed",   "", "", "response component (original) count")
header_rows$response_comp_cnt       <- c("", "label", "observed",   "", "", "response component count")
header_rows$sig_subm_id             <- c("", "label", "background", "", "", "sequential ID of signature within a publication (PMID) and for its response_type ")
header_rows$sig_row_id              <- c("", "label", "background", "", "", "row ID of signature within its submission data file(s)")
header_rows$row_key                 <- c("", "label", "background", "", "", "row key")

# header_rows and df2 should have the same column names (eventually)
df2 <- insub[first_data_row:nrow(insub),]

# unique observation ID (row number in original template) within this sheet
df2$sig_row_id  <- seq(from = first_data_row + 1, length.out = nrow(df2))

# The "pmid:" tag is not currently required, but check if anything was entered
df2$publication_reference_id <- sub("pmid:", "", df2$publication_reference_id)
s <- sapply(as.integer(df2$publication_reference_id), is.integer)
if (!all(s)) {
  stop("unexpected namespace tag")
}

# Infection templates contain more than one response_behavior_type, filter out all but current type
df2 <- df2[df2$response_behavior_type == response_behavior_type_var, ]

# remove signatures from unusable PMIDs
df2 <- df2[!(df2$publication_reference_id %in% exclude_pmid), ]
message("number of rows in df2: ", nrow(df2))

# Generate observation IDs based on the PMID field value and on original row number,
# to allow e.g. response_components to be summarized by observation later.
# These are generated before any row deletions so that can refer back to original template

# ID number of observation within its submission (based on PMID)
df2$sig_subm_id  <- cumcount(df2$publication_reference_id)
df2$row_key      <- paste(df2$publication_reference_id, df2$sig_subm_id, df2$sig_row_id, sep = "_")

#  This time, titles_and_dates_df will only have the PMIDs for the current exposure_type
titles_and_dates_df <- get_titles_and_dates(df2, RENEW_PMIDS, pmid_file, log_files, base_filename)
message("number of titles/dates: ", nrow(titles_and_dates_df))

# Most values are empty, not "Y" or "N".
# Currently only used for type GENE
df2$is_model[df2$is_model == ""] <- "N"

# fix capitalization (should fix in source google sheet)
df2$response_behavior_type <- tolower(df2$response_behavior_type)
df2$response_behavior      <- tolower(df2$response_behavior)

# check if  comment column has more than 255 character limit
w <- sapply(df2$comments, function(x) {nchar(x) > 255}, USE.NAMES	= FALSE)
if(any(w)) {
  df2$comments[w] <- stringr::str_trunc(df2$comments[w], 250)
  print(paste("some comments were too long for row(s), truncated:", paste(df2$row_key[w], collapse = ", ")))
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

codes <- strsplit(df2$exposure_material_id,";")  # returns a list if multiple values per row
# in Infection template, have allowed extra text in exposure_material_id column, need to handle
codes <- sub(" .*$", "", codes)
codes <- sub("ncbi_taxid:", "", codes)
s <- sapply(as.integer(codes), is.integer)
if (!all(s)) {
  print("unexpected namespace tag")
}
# this only works because only one entry per row
df2$exposure_material_id <- codes

############################################################################
##### Data substitution and splitting (1): response_component_original #####
############################################################################

# df2$response_component_original values become factors!
df2 <- cSplit(df2, "response_component_original", sep = ",", direction = "long")
df2 <- cSplit(df2, "response_component_original", sep = ";", direction = "long")

# using " /// " leaves extra slashes, regexp "[/]{3}" does not
# FIXME - getting warnings since updated to R 4.1:
# Warning message:
#   In type.convert.default(unlist(x, use.names = FALSE)) :
#   'as.is' should be specified by the caller; using TRUE
df2 <- cSplit(df2, "response_component_original", sep = "[/]{3}", direction = "long", fixed = FALSE)
# the curated data does in places have spaces as separators

# The cSplit() calls produce a data.table of data.frame rows.  Coerce back to just data.frame.
df2 <- as.data.frame(df2)

# remove factors from df2$response component so that for genes, can alter values.
# FIXME - any way to gain better control of this?
df2$response_component_original <- trimws(as.character(df2$response_component_original))

# create original gene symbols "signature" list after applying manual corrections

# Apply manual gene corrections
rvl <- manual_gene_corrections(df2$response_component_original,
                                 paste(reference_files, manual_gene_corrections_txt, sep = "/"))
genes <- rvl$genes

# Fix certain gene symbols containing "orf"
genes <- fix_orf_symbols(genes)

uids_list <- unique(df2$sig_row_id)
# get original count of rows for each sig_row_id.
uids_cnt  <- sapply(uids_list, function(x) {sum(df2$sig_row_id == x)})

for (i in 1:length(uids_list)) {
  # NOTE - df2$response_comp_orig_cnt is of type "character".  Be careful!
  df2$response_comp_orig_cnt[df2$sig_row_id == uids_list[i]] <- uids_cnt[i]
}

# Update original gene symbols to latest versions (NCBI/HGNC)

# starts with "genes" manually corrected list from above
rvl <- update_gene_symbols(genes, log_files, base_filename,
                             reference_files,
                             DOWNLOAD_NEW_HGNC)
genes_map <- rvl$genes_map
rm(rvl)

# Save the gene map for the unmatched symbols with their sig_row_id
# to add them back to complete signatures later
unmatched_symbols_map <- cbind(genes_map, sig_row_id = df2$sig_row_id)
unmatched_symbols_map <- unmatched_symbols_map[is.na(unmatched_symbols_map$Symbol), ]

# copy the fixed symbols back to the main data structure df2
df2$response_component <- genes_map$Symbol

# get rid of genes that had no valid symbol.
df2 <- df2[!is.na(df2$response_component), ]

# Remove exact duplicates (implies duplicate copies of response component in a signature)
# This can happen when e.g. two or more original probesets get map to one symbol.
df2 <- unique(df2)

# get count of rows for each sig_row_id.
uids_list <- unique(df2$sig_row_id)
uids_cnt  <- sapply(uids_list, function(x) {sum(df2$sig_row_id == x)})

for (i in 1:length(uids_list)) {
  df2$response_comp_cnt[df2$sig_row_id == uids_list[i]] <- uids_cnt[i]
}

#################################################
##### Data splitting (2): exposure material #####
#################################################

df2 <- cSplit(df2, "exposure_material_id", sep = ";", direction = "long")

if(any(df2$exposure_material_id == "")) {
  # Do not expect this
  stop("missing exposure_material_id")
}
# In case we forget to stop...
df2 <- df2[df2$exposure_material_id != "", ]

# The cSplit() calls produce a data.table of data.frame rows.  Coerce back to just data.frame.
df2 <- as.data.frame(df2, stringsAsFactors = FALSE)

###############################################################
##### Data splitting (4): tissue_type_term_id             #####
###############################################################

##### split ######
df2 <- cSplit(df2, "tissue_type_term_id", sep = ";", direction = "long")

# class(df2)  # "data.table" "data.frame"
df2 <- as.data.frame(df2)

df2$tissue_type_term_id <- sub(" .*$", "", df2$tissue_type_term_id)

###############################################################
##### Data splitting (5): comparison                      #####
###############################################################

##### split ######
df2 <- cSplit(df2, "comparison", sep = ";", direction = "long")

# class(df2)  # "data.table" "data.frame"
df2 <- as.data.frame(df2)
df2$comparison <- trimws(df2$comparison)

####################################################################
#### Get counts of each response component for e.g. word cloud #####
####################################################################

response_df <- unique(df2[ , c("response_component", "sig_row_id")])
response_df <- as.data.frame(table(response_df$response_component))
colnames(response_df) <- c("response_component", "count")
response_df <- response_df[order(response_df$count, decreasing = TRUE), ]

# Write out counts in tab-delimited format
# Use .txt extension rather than .tsv to thwart Excel auto-change of symbols
write.table(response_df,
            file = logfile_path(log_files, base_filename, "response_component_counts.txt"),
            sep = "\t", row.names = FALSE)

#############################################################
#### Write out full denormalized data                    ####
#############################################################
del_cols <- c("submission_name", "submission_date", "template_name", "short_comment")
df2tmp <- df2[!colnames(df2) %in% del_cols]

df2tmp <- df2tmp[-1]

filename <- logfile_path(standardized_curations, base_filename, "standardized_denormalized.tsv")
write.table(df2tmp,
            file = filename,
            sep = "\t", row.names = FALSE, col.names = TRUE)
gzip(filename, destname=paste0(filename, ".gz"), overwrite=TRUE, remove=TRUE)

#############################################################
#### Recreate original spreadsheet with all corrections #####
#############################################################

# Collect multi-valued columns by unique observation ID (each row in original spreadsheet)
uniq_sig_row_ids          <- unique(df2$sig_row_id)
resp_components_cnt_df    <- data.frame()
resp_components_annotated <- vector("list", length(uniq_sig_row_ids))
resp_components_collected <- vector("list", length(uniq_sig_row_ids))
resp_components_full_sig  <- vector("list", length(uniq_sig_row_ids))   # used for cell types only

for (i in 1:length(uniq_sig_row_ids)) {
  df2tmp <- df2[df2$sig_row_id == uniq_sig_row_ids[i], ]
  # Recreate a full signature in one row
  base_row <- df2tmp[1, ] # get first row for this uniqID

  response_rowname     <- paste(base_row$publication_reference_id, base_row$sig_subm_id, uniq_sig_row_ids[i], sep = "_")
  response_description <- paste("PMID", base_row$publication_reference_id, response_behavior_type_var, base_row$sig_subm_id, sep = " ")
  
  # Use the full original set of response components rather than just those
  # for which a valid symbol was found.
  base_row$response_component_original <- paste(unique(df2tmp$response_component_original), collapse = "; ")

  # The "updated" (gene or cell types after fixes) response_components
  resp_components_collected[[i]] <- unique(df2tmp$response_component)  # genes or base cell types

  base_row$exposure_material_id  <- paste(unique(df2tmp$exposure_material_id), collapse = "; ")
  base_row$tissue_type_term_id   <- paste(unique(df2tmp$tissue_type_term_id), collapse = "; ")
  
  
    base_row$response_component <- paste(unique(df2tmp$response_component), collapse = "; ")
    resp_components_annotated[[i]] <- c(response_rowname, response_description, unique(df2tmp$response_component))

    w <- which(titles_and_dates_df$pmid == base_row$publication_reference_id)
    if (length(w) != 1) {
      stop(paste("unexpected PMID result: uniqID = ", uniq_sig_row_ids[i], ", pmid = ", base_row$publication_reference_id))
    }
    titles_and_dates_row <- titles_and_dates_df[w, ]

  tmp <- data.frame(rowname = response_rowname,
                    pmid = base_row$publication_reference_id,
                    sig_subm_id = base_row$sig_subm_id,
                    sig_row_id = uniq_sig_row_ids[i],
                    count = length(resp_components_collected[[i]]))
  resp_components_cnt_df <- rbind(resp_components_cnt_df, tmp)
}

names(resp_components_collected) <- uniq_sig_row_ids  # name not actually used again?
names(resp_components_full_sig)  <- uniq_sig_row_ids  # name not actually used again?
names(resp_components_annotated) <- uniq_sig_row_ids  # name is used for this one.

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

df2$exposure_material <- NULL
# the real purpose of this script. the function name is misleading: this writes templates and all other release files.
write_submission_template(df2, header_rows, "../data/release_files", csv_submission_files, "hipc_inf_gene", titles_and_dates_df,
                            resp_components_collected, unmatched_symbols_map,
                            "GENE", "INFECTION", "Gene expression response to infection")

write.csv(resp_components_cnt_df,
          file = logfile_path(log_files, base_filename, "response_component_counts_by_row.csv"),
          row.names = FALSE)

# totals by PMID
rc_cnt_by_pmid <- lapply(unique(df2$publication_reference_id), function(pmid) {
  vals <- resp_components_cnt_df[resp_components_cnt_df$pmid == pmid, "count"]
  return(list(pmid = pmid, sum = sum(vals)))
})
rc_cnt_by_pmid <- rbindlist(rc_cnt_by_pmid)
# Note that response components can appear in more than one signature per PMID.  These are not unique counts.
write.csv(rc_cnt_by_pmid,
          file = logfile_path(log_files, base_filename, "response_component_non-unique_count_by_PMID.csv"),
          row.names = FALSE)

# Unique response component count by PMID
# ordered by count, but could as well be ordered by PMID
rc_cnt_by_pmid_uniq <- resp_comps_by_pmid(resp_components_annotated, unique(df2$publication_reference_id))
rc_cnt_by_pmid_uniq <- rc_cnt_by_pmid_uniq[order(rc_cnt_by_pmid_uniq$count, decreasing = TRUE), ]

write.csv(rc_cnt_by_pmid_uniq,
          file = logfile_path(log_files, base_filename, "response_component_unique_count_by_PMID.csv"),
          row.names = FALSE)
