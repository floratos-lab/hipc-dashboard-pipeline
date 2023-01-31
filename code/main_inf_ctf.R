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
source("pmid_to_title_easy.R")
source("write_submissions.R")

# Assume executing from the ./code directory
source_curations       <- "../data/source_curations"
reference_files        <- "../data/reference_files"

log_files              <- "../logfiles"

ctf_fixes_tsv               <- "cell_type_frequency-response_components_mapping.txt"

exclude_pmid <- "33361761"  # PMID(s) that are not suitable for processing

##### Set up file and template name components #####
first_data_row <- 9  # not counting the first row, which becomes a column header.
sheet_file     <- "COVID-19 curation template - example curation.tsv"
sheet_file2    <- "Odak_2020_pmid_32650275 - covid19.tsv"
sheet_file3    <- "hipc_infection_covid_v2 - multiple_types.tsv"

base_filename  <- "inf_cell_type"

# used to filter sheet rows
response_behavior_type_var <- "cell-type frequency"

COLUMN_NAMES = c("X", "submission_name", "submission_date", "template_name", "response_component",
  "curation_date", "cohort", "age_min", "age_max", "age_units", "number_subjects", "tissue_type", 
  "tissue_type_term_id", "method", "response_component_original", "is_model", "response_behavior_type", 
  "response_behavior", "comparison", "baseline_time_event", "time_point", "time_point_units", 
  "exposure_material_id", "exposure_process", "disease_name", "disease_stage", "publication_reference_id", 
  "publication_date", "publication_reference_url", "signature_source", "comments", 
  "response_comp_orig_cnt", "response_comp_cnt", "sig_subm_id", "sig_row_id", "row_key")

insub <- read.delim(file =  paste(source_curations, sheet_file, sep = "/"),
                    strip.white = TRUE,
                    stringsAsFactors = FALSE)
insub <- insub[colnames(insub) %in% COLUMN_NAMES]

# read two other source curation files
insub2 <- read.delim(file =  paste(source_curations, sheet_file2, sep = "/"),
                       strip.white = TRUE,
                       stringsAsFactors = FALSE)
insub2 <- insub2[colnames(insub2) %in% COLUMN_NAMES]

insub3 <- read.delim(file =  paste(source_curations, sheet_file3, sep = "/"),
                         strip.white = TRUE,
                         stringsAsFactors = FALSE,
                         quote = "")
insub3 <- insub3[colnames(insub3) %in% COLUMN_NAMES]

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
insub <- data.frame(insub[, 1],
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
message("number of rows in the ctf mapping file ", nrow(ctf_fixes))

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

failed_matches <- df2$response_component_original[is.na(ctf_match)]
message("unique response_component not in the ctf mapping file ", length(unique(failed_matches)))
successful_matches <- df2$response_component_original[!is.na(ctf_match)]
message("number of rows with response_component in the ctd mapping file ", length(successful_matches))

# combine proterm(s) and extra together.  All same length as df2.
proterm <- ctf_fixes$PRO_term_w_intensity[ctf_match]
extra <- ctf_fixes$extra[ctf_match]
df2$response_component <- ctf_fixes$CL_term[ctf_match]
df2$response_component_id   <- ctf_fixes$CL_ID[ctf_match]
df2$proterm_and_extra  <- ifelse(proterm!="" & extra!="", paste0("& ", proterm, ", ", extra), 
    ifelse(proterm!="", paste("&", proterm), "")
)
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

write_submission_template(df2, header_rows, "../data/release_files", NULL, "hipc_inf_ctf", titles_and_dates_df,
                            resp_components_full_sig, unmatched_symbols_map = NULL,
                            "CELLTYPE_FREQUENCY", "INFECTION", "Immune cell-type frequency response to infection")
