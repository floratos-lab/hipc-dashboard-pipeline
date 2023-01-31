# Description:
#  Starting with curated HIPC data, prepare Dashboard submission templates.
#  Type: Cell type frequency.
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
#   * "cell_type_frequency-response_components_mapping.txt" - corrections to exposure_material
#        terms for cell type frequency sheet.  These corrections will in the end be added
#        directly into the HIPC Dashboard spreadsheet.
#
# Output files:
#   * Submission file for each signature
#   * Files containing the list of response_components for each signature.
#       One file has only the HGNC symbols, the other has all original symbols.
#   * CV-per-template file for all signatures of a type (gene, cell-type)

library(splitstackshape)  # for cSplit()
library(uniqtag)          # for cumcount()

source("hipc_utils.R")
source("vaccines_pathogens.R")
source("pmid_to_title_easy.R")
source("write_submissions.R")

# Assume executing from the ./code directory
source_curations       <- "../data/source_curations"
reference_files        <- "../data/reference_files"

log_files              <- "../logfiles"

vaccine_tsv                 <- "vaccine_years.txt"
ctf_fixes_tsv               <- "cell_type_frequency-response_components_mapping.txt"

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
exclude_pmid <- "33361761"  # PMID(s) that are not suitable for processing

##### Set up file and template name components #####
# change sheet name spaces to underscores
first_data_row <- 7  # not counting the first row, which becomes a column header.

sheet_file     <- "hipc_vaccine - cell_type_frequency.tsv"
base_filename  <- "vac_cell_type"

# used to filter sheet rows
response_behavior_type_var <- "cell-type frequency"

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
                    response_component_id  = "",
                    proterm_and_extra = "",
                    pro_ontology_id = "",
                    fully_qualified_response_component = "",
                    response_comp_orig_cnt = "",
                    response_comp_cnt = "",
                    sig_subm_id = "",
                    sig_row_id = "",
                    row_key = "",
                    stringsAsFactors = FALSE)

  # cell-type specific headers
  # Need to set for INFECTION type, will just overwrite existing for VACCINE
header_rows <- insub[1:6,]
header_rows$response_component      <- c("", "label", "observed", "", "", "response component name")
header_rows$response_component_id   <- c("cell_subset", "", "cell_biomarker", "", "", "response component")
header_rows$proterm_and_extra       <- c("", "label", "observed", "", "", "response component protein ontology and extra terms")
header_rows$pro_ontology_id         <- c("", "label", "observed", "", "", "response component marker protein ontology ID")
header_rows$fully_qualified_response_component <- c("", "label", "observed", "", "", "response component (fully qualified)")
header_rows$response_component_original <- c("", "label", "observed", "", "", "response component (original curated cell type)")
header_rows$submission_name         <- c("", "label", "background", "", "", "submission name")
header_rows$submission_date         <- c("", "label", "background", "", "", "submission_date")
header_rows$template_name           <- c("", "label", "background", "", "", "template_name")
header_rows$response_comp_orig_cnt  <- c("", "label", "observed",   "", "", "response component (original) count")
header_rows$response_comp_cnt       <- c("", "label", "observed",   "", "", "response component count")
header_rows$sig_subm_id             <- c("", "label", "background", "", "", "sequential ID of signature within a publication (PMID) and for its response_type ")
header_rows$sig_row_id              <- c("", "label", "background", "", "", "row ID of signature within its submission data file(s)")
header_rows$row_key                 <- c("", "label", "background", "", "", "row key")

ctf_fixes <- read.delim(file = paste(reference_files, ctf_fixes_tsv, sep = "/"),
                          stringsAsFactors = FALSE)
message("dim(ctf_fixes) = ", dim(ctf_fixes))

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
text  <- strsplit(df2$exposure_material, ";") # returns a list
m <- mapply(function(x, y) {length(x) == length(y)}, codes, text)
if (!all(m)) {
  warning(paste("exposure code vs text length mismatch for rows rowkeys:", paste(df2$row_key[!m], collapse = ", ")))
  warning(paste("removing non-matching rows"))
  df2 <- df2[m, ]
  codes <- codes[m]
  text <- text[m]
}

############################################################################
##### Data substitution and splitting (1): response_component_original #####
############################################################################

df2 <- cSplit(df2, "response_component_original", sep = ";", direction = "long")


# The cSplit() calls produce a data.table of data.frame rows.  Coerce back to just data.frame.
df2 <- as.data.frame(df2)

# remove factors from df2$response component so that for genes, can alter values.
# FIXME - any way to gain better control of this?
df2$response_component_original <- trimws(as.character(df2$response_component_original))

uids_list <- unique(df2$sig_row_id)
# get original count of rows for each sig_row_id.
uids_cnt  <- sapply(uids_list, function(x) {sum(df2$sig_row_id == x)})

for (i in 1:length(uids_list)) {
  # NOTE - df2$response_comp_orig_cnt is of type "character".  Be careful!
  df2$response_comp_orig_cnt[df2$sig_row_id == uids_list[i]] <- uids_cnt[i]
}

# Update original gene symbols to latest versions (NCBI/HGNC)

# ctf_match will have NA values where no match found
ctf_match <- match(df2$response_component_original, ctf_fixes$original_annotation)

# combine proterm(s) and extra together.  All same length as df2.
proterm <- ctf_fixes$PRO_term_w_intensity[ctf_match]
extra <- ctf_fixes$extra[ctf_match]
# vector initialized with ""
proterm_and_extra <- vector(mode = "character", length = length(proterm))

for (i in 1:length(proterm)) {
    if(is.na(ctf_match[i])) {
      next
    } 
    # insert comma between non-empty terms
    a <- c(proterm[i], extra[i])
    a <- a[a!=""]                   # returns char(0) if both empty
    a <- paste(a, collapse = ", ")  # returns "" on empty
    proterm_and_extra[i] <- ifelse(a != "", paste("&", a), "")
  }

  df2$response_component <- ctf_fixes$CL_term[ctf_match]
  df2$response_component_id   <- ctf_fixes$CL_ID[ctf_match]
  df2$proterm_and_extra  <- proterm_and_extra
  df2$pro_ontology_id    <- ctf_fixes$PRO_ID[ctf_match]
  df2$fully_qualified_response_component <- ""

if(nrow(df2) != length(ctf_match)) {
  stop("nrow(df2) != length(ctf_match)")
}
# Remove any unmatched cell types (should have already been caught)
df2 <- df2[!is.na(ctf_match), ]

concatenated_cell_types <- trimws(paste(df2$response_component, df2$proterm_and_extra, sep = " "))
df2$fully_qualified_response_component <- concatenated_cell_types

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
resp_components_full_sig  <- vector("list", length(uniq_sig_row_ids))   # used for cell types only

for (i in 1:length(uniq_sig_row_ids)) {
  df2tmp <- df2[df2$sig_row_id == uniq_sig_row_ids[i], ]

  # FIXME - should not need unique here, still just one component per row
  full_sig <- unique(df2tmp$fully_qualified_response_component)
  # The pro_ontology_id values are already separated by semicolons, so change to commas
  # before potentially joining two lists of pro-terms.  
  df2tmp$pro_ontology_id <- sapply(df2tmp$pro_ontology_id, function(x) {gsub(";", ",", x)})
    
  resp_components_full_sig[[i]]  <- full_sig
}

names(resp_components_full_sig)  <- uniq_sig_row_ids  # name not actually used again?

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

write_submission_template(df2, header_rows, "../data/release_files", NULL, "hipc_vac_ctf", titles_and_dates_df,
                            resp_components_full_sig, unmatched_symbols_map = NULL,
                            "CELLTYPE_FREQUENCY", "VACCINE", "Immune cell-type frequency response to vaccine exposure")
