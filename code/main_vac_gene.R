# Description:
#  Starting with curated HIPC data, prepare Dashboard submission templates.
#  Type: Gene
#  * Perform gene symbol updates to conform to HGNC.
#    * Fix up temporary term problems pending e.g. ontology updates
#    * When exposure material is an influenza vaccine,
#        substitute the actual virus components for the listed vaccine year
#        into the target pathogen field.
#    * Retrieve publication title and abstract etc. based on PMID
#    * Create individual submission files by publication
#    * Create CV-per-template file for set of submissions
#
# Input files (must be in current directory):
#   * "vaccine_years.txt" - maps vaccine season year to vaccine viral components
#   * "manual_gene_symbol_corrections.txt" - maps invalid symbols  to known valid sybmols
#
# Output files:
#   * Submission file for each signature
#   * Files containing the list of response_components for each signature.
#       One file has only the HGNC symbols, the other has all original symbols.
#   * CV-per-template file for all signatures of a type (gene, cell-type)

library(splitstackshape)  # for cSplit()
library(uniqtag)          # for cumcount()

source("hipc_utils.R")
source("gene_routines.R")
source("vaccines_pathogens.R")
source("pmid_to_title_easy.R")
source("write_submissions.R")

# Assume executing from the ./code directory
source_curations       <- "../data/source_curations"
reference_files        <- "../data/reference_files"
log_files              <- "../logfiles"

vaccine_tsv                 <- "vaccine_years.txt"
manual_gene_corrections_txt <- "manual_gene_symbol_corrections.txt"

# In the observation summary, do not display pathogens if the vaccine already
# uses the pathogen in its name:
# VO_0004810: 2011-2012 trivalent inactivated vaccine (A/California/7/09 (H1N1,),
#     A/Perth /16/2009 (H3N2), and B/Brisbane/60/2008).
# VO_0004899: 2012-2013 seasonal trivalent inactivated influenza vaccine
#     (A/California/7/2009 (H1N1), A/Victoria/361/2011 (H3N2),
#     and B/Wisconsin/1/2010)
# VO_0004903: Inactivated monovalent influenza A/H5N1
#     (3.75 mcg hemagglutinin [HA] A/Indonesia/05/2005) split-virus (SV) vaccine (Sanofi)
vaccine_VO_has_pathogens <- c("VO_0004810", "VO_0004899", "VO_0004903") # global variable used in write_submissions, not locally
exclude_pmid <- "33361761"  # PMID(s) that are not suitable for processing

##### Set up file and template name components #####
first_data_row <- 7  # not counting the first row, which becomes a column header.

sheet_file     <- "hipc_vaccine - gene_expression.tsv"
base_filename  <- "vac_gene_expression"

# used to filter sheet rows
response_behavior_type_var <- "gene expression"

vaccines_by_year <- read.delim(file = paste(reference_files, vaccine_tsv, sep = "/"),
                                 stringsAsFactors = FALSE)

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

del_cols <- c("method")  # not filled in
insub <- insub[!(colnames(insub) %in% del_cols)]

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
# Need to set for INFECTION type, will just overwrite existing for VACCINE
header_rows = insub[1:6,]
header_rows$response_component <- c("gene", "", "gene_biomarker", "", "", "response component")
header_rows$response_component_original <- c("", "label", "observed", "", "", "response component (original gene symbol)")
header_rows$submission_name    <- c("", "label", "background", "", "", "submission name")
header_rows$submission_date    <- c("", "label", "background", "", "", "submission_date")
header_rows$template_name      <- c("", "label", "background", "", "", "template_name")
header_rows$response_comp_orig_cnt  <- c("", "label", "observed",   "", "", "response component (original) count")
header_rows$response_comp_cnt  <- c("", "label", "observed",   "", "", "response component count")
header_rows$sig_subm_id        <- c("", "label", "background", "", "", "sequential ID of signature within a publication (PMID) and for its response_type ")
header_rows$sig_row_id         <- c("", "label", "background", "", "", "row ID of signature within its submission data file(s)")
header_rows$row_key            <- c("", "label", "background", "", "", "row key")

df2 <- insub[first_data_row:nrow(insub),]

# unique observation ID (row number in original template) within this sheet
df2$sig_row_id  <- seq(from = first_data_row + 1, length.out = nrow(df2))

# The "pmid:" tag is not currently required, but check if anything was entered
df2$publication_reference_id <- sub("pmid:", "", df2$publication_reference_id)
s <- sapply(as.integer(df2$publication_reference_id), is.integer)
if (!all(s)) {
  stop("unexpected namespace tag")
}

df2 <- df2[df2$publication_reference_id != "", ]

# remove signatures from unusable PMIDs
df2 <- df2[!(df2$publication_reference_id %in% exclude_pmid), ]
message("number of rows in df2: ", nrow(df2))

# Generate observation IDs based on the PMID field value and on original row number,
# to allow e.g. response_components to be summarized by observation later.
# These are generated before any row deletions so that can refer back to original template

# ID number of observation within its submission (based on PMID)
df2$sig_subm_id  <- cumcount(df2$publication_reference_id)
df2$row_key      <- paste(df2$publication_reference_id, df2$sig_subm_id, df2$sig_row_id, sep = "_")

# get titles_and_dates_df
load(file = paste0(reference_files, "/", base_filename, "-titles_and_dates_df.RData"))
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

##########################################################
##### Other global changes to template before splits #####
##########################################################

## Create a map of VO codes to text equivalents
## Both columns must have the same number of values per row
codes <- strsplit(df2$exposure_material_id,";")  # returns a list if multiple values per row
text  <- strsplit(df2$exposure_material, ";") # returns a list
m <- mapply(function(x, y) {length(x) == length(y)}, codes, text)
if (!all(m)) {
  warning(paste("exposure code vs text length mismatch for rows rowkeys:", paste(df2$row_key[!m], collapse = ", ")))
  warning(paste("removing non-matching rows"))
  df2 <- df2[m, ]
}

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
rvl <- update_gene_symbols(genes, log_files, base_filename, reference_files)
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
##### Data splitting (3): target pathogens (VACCINE ONLY) #####
###############################################################

# Substitute in actual virus components for influenza pathogens only.
# Note - this is done before the real splitting below because the
# lookup_vaccine() function needs all codes at once for an observation.
# lookup_vaccine() depends on global variable 'vaccines_by_year'!!
for (i in 1:nrow(df2)) {
  if(grepl("influ:", df2$target_pathogen_taxonid[i])) {
    codes <- unlist(strsplit(df2$target_pathogen_taxonid[i],";"))
    codes <- sub("influ:", "", codes)
    codes <- trimws(codes)
    codes <- paste(codes, collapse = ", ")
    df2$target_pathogen[i]         <- lookup_vaccine(codes, type = "name")
    df2$target_pathogen_taxonid[i] <- lookup_vaccine(codes, type = "taxid")
  }
}

# Substitute in actual pathogens for measles + DPT vaccination
# Note - the specific pathogen strains below are not specified by the VO codes used.
# These are the only case in Gene and Cell-type where there are more than a single pathogen,
# other than for influenza, so have to rewrite the pathogen entries to match split vaccine VO codes
w <- which(df2$exposure_material_id == "VO_0000654")
df2[w, "target_pathogen"] <- "Measles virus strain Edmonston-Zagreb"
df2[w, "target_pathogen_taxonid"] <- "70149"
w <- which(df2$exposure_material_id == "VO_0000738")
df2[w, "target_pathogen"] <- "Corynebacterium diphtheriae; Clostridium tetani; Human poliovirus 1; Human poliovirus 2; Human poliovirus 3"
df2[w, "target_pathogen_taxonid"] <- "1717; 1513; 12080; 12083; 12086"

##### split ######
df2 <- cSplit(df2, "target_pathogen_taxonid", sep = ";", direction = "long")
df2 <- as.data.frame(df2)

# in Infection template, have allowed "(name)" exposure_material_id column.  Remove.
codes <- sub(" .*$", "", df2$target_pathogen_taxonid)
codes <- sub("ncbi_taxid:", "", codes)
s <- sapply(as.integer(codes), is.integer)
if (!all(s)) {
  print("target_pathogen_taxonid: unexpected namespace tag")
}
df2$target_pathogen_taxonid <- codes

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

#############################################################
#### Recreate original spreadsheet with all corrections #####
#############################################################

# Collect multi-valued columns by unique observation ID (each row in original spreadsheet)
uniq_sig_row_ids          <- unique(df2$sig_row_id)
resp_components_collected <- vector("list", length(uniq_sig_row_ids))

for (i in 1:length(uniq_sig_row_ids)) {
  df2tmp <- df2[df2$sig_row_id == uniq_sig_row_ids[i], ]
  # Recreate a full signature in one row
  base_row <- df2tmp[1, ] # get first row for this uniqID

  # Use the full original set of response components rather than just those
  # for which a valid symbol was found.
  base_row$response_component_original <- paste(unique(df2tmp$response_component_original), collapse = "; ")

  # The "updated" (gene or cell types after fixes) response_components
  resp_components_collected[[i]] <- unique(df2tmp$response_component)  # genes or base cell types

  base_row$exposure_material_id  <- paste(unique(df2tmp$exposure_material_id), collapse = "; ")
  base_row$tissue_type_term_id   <- paste(unique(df2tmp$tissue_type_term_id), collapse = "; ")
  base_row$response_component <- paste(unique(df2tmp$response_component), collapse = "; ")

  w <- which(titles_and_dates_df$pmid == base_row$publication_reference_id)
  if (length(w) != 1) {
    stop(paste("unexpected PMID result: uniqID = ", uniq_sig_row_ids[i], ", pmid = ", base_row$publication_reference_id))
  }
}

names(resp_components_collected) <- uniq_sig_row_ids  # name not actually used again?

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
del_cols <- c("exposure_material", "target_pathogen", "short_comment", "process_note")

df2 <- df2[!colnames(df2) %in% del_cols]
header_rows <- header_rows[!colnames(header_rows) %in% del_cols]

write_submission_template(df2, header_rows, "../data/release_files", NULL, "hipc_vac_gene", titles_and_dates_df,
                            resp_components_collected, unmatched_symbols_map,
                            "GENE", "VACCINE", "Gene expression response to vaccine exposure")
