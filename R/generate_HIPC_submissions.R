# File name: generate_HIPC_submissions.R
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
#
# BiocManager::install(version = "3.14")
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
# BiocManager::install("stringr")
# BiocManager::install("pracma")

library(splitstackshape)  # for cSplit()
library(uniqtag)          # for cumcount()
# Note that rbindlist() tends to return a data.table with a data.frame, which causes
# no end of problems if not cast back to just a data.frame.
library(data.table)       # for rbindlist()
library(ggplot2)
library(stringi)
library(stringr)


source("hipc_utils.R")
source("gene_routines.R")
source("vaccines_pathogens.R")
source("pmid_to_title_easy.R")
source("write_submissions.R")
source("msigdb_submission_utils.R")
source("find_unique.R")

# Available response_type values are "GENE", "CELLTYPE_FREQUENCY"
response_type <- "GENE"
# Available exposure_type values are "VACCINE", "INFECTION" (covid-19)
exposure_type <- "VACCINE"

# For the moment, assume executing interactively from the ./R directory
source_data_dir <- "../source_data"
submission_dir  <- "../submissions"
logdir          <- "../logfiles"

vaccine_tsv     <- "vaccine_years.txt"
ctf_fixes_tsv   <- "cell_type_frequency-response_components_mapping.txt"
manual_gene_corrections_tsv <- "manual_gene_symbol_corrections.txt"

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

if (response_type == "GENE" && DOWNLOAD_NEW_NCBI_GENES) {
  if(!update_ncbi_homo_sapiens_gene_info(source_data_dir)) {
    print("update of NCBI gene info file failed")
  }
}

options(stringsAsFactors = FALSE)  # unfortunately doesn't help with cSplit output
summary_df <- data.frame()  # initialize summary log

##### Set up file and template name components #####
# change sheet name spaces to underscores
if (response_type == "GENE" && exposure_type == "VACCINE") {
  sheet_file     <- "hipc_vaccine - gene_expression.tsv"
  base_filename  <- "vac_gene_expression"
  template_name  <- "hipc_vac_gene"
  project        <- "Gene expression response to vaccine exposure"
} else if (response_type == "GENE" && exposure_type == "INFECTION") {
  sheet_file     <- "COVID-19 curation template - example curation.tsv"
  sheet_file2    <- "Odak_2020_reappearance_of_effector_pmid_32650275 - covid19.tsv"
  sheet_file3    <- "hipc_infection_covid_v2 - multiple_types.tsv"
  base_filename  <- "inf_gene_expression"
  template_name  <- "hipc_inf_gene"
  project        <- "Gene expression response to infection"
} else if (response_type == "CELLTYPE_FREQUENCY" && exposure_type == "VACCINE") {
  sheet_file     <- "hipc_vaccine - cell_type_frequency.tsv"
  base_filename  <- "vac_cell_type"
  template_name  <- "hipc_vac_ctf"
  cell_mapping_sheet_name   <- "HIPC_Dashboard-Cell_type_Freque"
  project        <- "Immune cell-type frequency response to vaccine exposure"
} else if (response_type == "CELLTYPE_FREQUENCY" && exposure_type == "INFECTION") {
  sheet_file     <- "COVID-19 curation template - example curation.tsv"
  sheet_file2    <- "Odak_2020-Julia_Davis-Porada-covid19.tsv"
  sheet_file3    <- "hipc_infection_covid_v2 - multiple_types.tsv"
  base_filename  <- "inf_cell_type"
  template_name  <- "hipc_inf_ctf"
  cell_mapping_sheet_name   <- "HIPC_Dashboard-Cell_type_Freque"
  project        <- "Immune cell-type frequency response to infection"
} else {
  stop("unknown sheet type")
}
# for joint summary, list all possible values of base_filename from above
if (exposure_type == "VACCINE") {
  all_response_types <- c("vac_gene_expression", "vac_cell_type")
} else if (exposure_type == "INFECTION") {
  all_response_types <- c("inf_gene_expression", "inf_cell_type")
}
pmid_file <- paste(source_data_dir,
                   paste(base_filename, "titles_and_dates_df.RData", sep = "-"),
                   sep = "/")

if (exposure_type == "VACCINE") {
  vaccines_by_year <- read.delim(file = paste(source_data_dir, vaccine_tsv, sep = "/"),
                                 stringsAsFactors = FALSE)
}

insub <- read.delim(file =  paste(source_data_dir, sheet_file, sep = "/"),
                    strip.white = TRUE,
                    stringsAsFactors = FALSE)
nrow(insub)

# Get rid of unused columns (empty in curated data) or those not meant to appear in the Dashboard.
# No error if named column does not exist in a particular sheet.
# Note - columns "short_comment" and "process_note" are removed later after not needed anymore
del_cols_common <- c("submission_name", "template_name", # should not appear anymore
              "spot_check", "Spot.check",
              "second_spot_check",
              "Curator",
              "curator_comments")
insub <- insub[!(colnames(insub) %in% del_cols_common)]

if (exposure_type == "VACCINE") {
  del_cols <- c("method")  # not filled in
  insub <- insub[!(colnames(insub) %in% del_cols)]

} else if(exposure_type == "INFECTION") {
  # read additional file
  # FIXME - a more general solution will be needed for handling multiple files
  insub2 <- read.delim(file =  paste(source_data_dir, sheet_file2, sep = "/"),
                      strip.white = TRUE,
                      stringsAsFactors = FALSE)
  nrow(insub2)
  insub2 <- insub2[!(colnames(insub2) %in% del_cols_common)]

#  insub3 <- read.delim(file =  paste(source_data_dir, sheet_file3, sep = "/"),
#                       strip.white = TRUE,
#                       stringsAsFactors = FALSE,
#                       quote = "")
#  nrow(insub3)
#  insub3 <- insub3[!(colnames(insub3) %in% del_cols_common)]

  del_cols <- c("route",
                "addntl_time_point_units",
                "signature_source_url")  # discontinued in latest template
  colnames(insub)[(colnames(insub) %in% del_cols)]

  colnames(insub2)[(colnames(insub2) %in% del_cols)]
  insub2 <- insub2[!(colnames(insub2) %in% del_cols)]

#  colnames(insub3)[(colnames(insub3) %in% del_cols)]
#  insub3 <- insub3[!(colnames(insub3) %in% del_cols)]

#  colnames(insub3)[!(colnames(insub3) %in% (colnames(insub2)))]
#  colnames(insub2)[!(colnames(insub2) %in% (colnames(insub3)))]

  all(colnames(insub) == colnames(insub2))
#  all(colnames(insub2) == colnames(insub3))

  # stitch together the data sections
  nrow(insub)
  insub <- rbind(insub, insub2[9:nrow(insub2),])
  nrow(insub)
#  insub <- rbind(insub, insub3[9:nrow(insub3),])
#  nrow(insub)

}

uniq_pmid <- unique(insub$publication_reference_id[9:nrow(insub)])
uniq_pmid <- sub("pmid:", "", uniq_pmid)
uniq_pmid <- uniq_pmid[uniq_pmid != ""]


### Make any changes to headers right at the beginning,
### before separating headers and data
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
if (response_type == "GENE") {
  # gene-specific headers
  # Need to set for INFECTION type, will just overwrite existing for VACCINE
  insub$response_component[1:6] <- c("gene", "", "gene_biomarker", "", "", "response component")
  insub$response_component_original[1:6] <-
    c("", "label", "observed", "", "", "response component (original gene symbol)")

} else if (response_type == "CELLTYPE_FREQUENCY") {
  insub$response_component_id  <- ""
  insub$proterm_and_extra <- ""
  insub$pro_ontology_id   <- ""
  insub$fully_qualified_response_component <- ""

  # cell-type specific headers
  # Need to set for INFECTION type, will just overwrite existing for VACCINE
  insub$response_component[1:6] <- c("", "label", "observed", "", "", "response component name")
  insub$response_component_id[1:6]  <-
    c("cell_subset", "", "cell_biomarker", "", "", "response component")
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

if (response_type == "CELLTYPE_FREQUENCY") {
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
if (exposure_type == "VACCINE") {  # the original curations, assumes no more will be added using updated template
  df2 <- insub[7:nrow(insub),]
} else {  # INFECTION, the new curation template has two extra rows of instructions
  df2 <- insub[9:nrow(insub),]
}

# Special handling for INFECTION templates
if(exposure_type == "INFECTION") {
  if (response_type == "CELLTYPE_FREQUENCY") {
    w <- df2$response_behavior_type == "cell-type frequency"
  } else if (response_type == "GENE") {
    w <- df2$response_behavior_type == "gene expression"
  } else {
    stop("unexpected sheet type, not yet implemented")
  }
  nrow(df2)
  df2 <- df2[w, ]
  nrow(df2)
}

# Depending on how the text file is created, the sheet may have blank rows at bottom.
# NOTE - the step above for INFECTION deals with this problem already
nrow(df2)
df2 <- df2[df2$publication_reference_id != "", ]
nrow(df2)
summary_df <- add_to_summary(summary_df, "Data rows in original sheet", nrow(df2))


# Generate observation IDs based on the PMID field value and on original row number,
# to allow e.g. response_components to be summarized by observation later.
# These are generated before any row deletions so that can refer back to original template
# The "pmid:" tag is not currently required, but check if anything was entered
df2$publication_reference_id <- sub("pmid:", "", df2$publication_reference_id)
s <- sapply(as.integer(df2$publication_reference_id), is.integer)
if (!all(s)) {
  stop("unexpected namespace tag")
}

# ID number of observation within its submission (based on PMID)
df2$subm_obs_id  <- cumcount(df2$publication_reference_id)
# unique observation ID (row number in original template) within this sheet
df2$uniq_obs_id  <- seq(from = 8, length.out = nrow(df2))
# save so can check later for removed observations
uniq_obs_id_orig <- df2$uniq_obs_id
df2$row_key      <- paste(df2$publication_reference_id, df2$subm_obs_id, df2$uniq_obs_id, sep = "_")

########################################################
##### Fix temporary problems with the curated data #####
########################################################

# Remove rows marked "skip"
# FIXME - not being used anymore in HIPC Dashboard source
# FIXME -should actually test for presence of this column
#        rather than it being template-dependent.
if(exposure_type == "VACCINE") {
  skip <- grepl("^skip", df2$process_note, ignore.case = TRUE)
  df2 <- df2[!skip, ]
}

# FIXME - Most values are empty, not "Y" or "N".
#         Empty is not necessarily the same as "N".
# Currently only used for type GENE
if(response_type == "GENE") {
  df2$is_model[df2$is_model == ""] <- "N"
}

# Write list of pathogens before any substitutions
if(exposure_type == "VACCINE") {
  write_unique_list(df2$target_pathogen, logdir, base_filename, "target_pathogens_before_fixes", do_split = TRUE)
}

# fix capitalization (should fix in source google sheet)
df2$response_behavior_type <- tolower(df2$response_behavior_type)
df2$response_behavior      <- tolower(df2$response_behavior)

# check if  comment column has more than 255 character limit
w <- sapply(df2$comments, function(x) {nchar(x) > 255}, USE.NAMES	= FALSE)
if(any(w)) {
  df2$comments[w] <- stringr::str_trunc(df2$comments[w], 250)
  print(paste("some comments were too long for row(s), truncated:", paste(df2$row_key[w], collapse = ", ")))
}

# check if comment column has non-UTF8 characters
# FIXME - this can't happen because going through text files.
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
codes <- strsplit(df2$exposure_material_id,";")  # returns a list if multiple values per row
# in Infection template, have allowed extra text in exposure_material_id column, need to handle
if(exposure_type == "INFECTION") {
  codes <- sub(" .*$", "", codes)
  codes <- sub("ncbi_taxid:", "", codes)
  s <- sapply(as.integer(codes), is.integer)
  if (!all(s)) {
    print("unexpected namespace tag")
  }
  # FIXME - this only works because only one entry per row
  df2$exposure_material_id <- codes
}

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

exposure_code_map <- unique(data.frame(codes, text, stringsAsFactors = FALSE))
exposure_code_map <- exposure_code_map[order(exposure_code_map$codes), ]

write.table(exposure_code_map,
            file = logfile_path(logdir, base_filename, "exposure_material_code_mapping.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE)

dupcodes <- exposure_code_map[duplicated(exposure_code_map$codes), ]$codes
w <- which(exposure_code_map$codes %in% dupcodes)
write.table(exposure_code_map[w, ],
            file = logfile_path(logdir, base_filename, "exposure_material_code_duplicates.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(unique(text),
            file = logfile_path(logdir, base_filename, "exposure_materials.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE)


############################################################################
##### Data substitution and splitting (1): response_component_original #####
############################################################################

if (response_type == "GENE") {
  # df2$response_component_original values become factors!
  df2 <- cSplit(df2, "response_component_original", sep = ",", direction = "long")
  df2 <- cSplit(df2, "response_component_original", sep = ";", direction = "long")
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
} else if (response_type == "CELLTYPE_FREQUENCY") {
  df2 <- cSplit(df2, "response_component_original", sep = ";", direction = "long")
}

# The cSplit() calls produce a data.table of data.frame rows.  Coerce back to just data.frame.
df2 <- as.data.frame(df2)

# remove factors from df2$response component so that for genes, can alter values.
# FIXME - any way to gain better control of this?
df2$response_component_original <- as.character(df2$response_component_original)


if (response_type == "GENE") {
  # check for lists of response components within one PMID that have a large overlap with one another.
  # This is intended to catch the case were one set was accidentally appended to another.
  # (It happens!)
  ft <- check_response_components_overlap(df2, unique(df2$publication_reference_id),
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
if (response_type == "GENE") {
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
} else if (response_type == "CELLTYPE_FREQUENCY") {

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
if (response_type == "GENE") {
  # starts with "genes" manually corrected list from above
  # make sure "genes" is still in synch with df2.
  length(genes) == nrow(df2)
  rvl <- update_gene_symbols(genes, logdir, base_filename,
                             source_data_dir,
                             DOWNLOAD_NEW_HGNC)
  genes_map <- rvl$genes_map
  summary_df <- rbind(summary_df, rvl$summary)
  rm(rvl)
  
  # Save PMID info for unmapped symbols
  no_valid_symbols_df <- df2[is.na(genes_map$Symbol), names(df2) %in%
                             c("response_component_original", "publication_reference_id", "subm_obs_id", "uniq_obs_id")]
  log_no_valid_symbol_vs_pmid(no_valid_symbols_df, base_filename)

  # w <- which(genes_map$Symbol != genes_map$alias)

  # Save the gene map for the unmatched symbols with their uniq_obs_id
  # to add them back to complete signatures later
  unmatched_symbols_map <- cbind(genes_map, uniq_obs_id = df2$uniq_obs_id)
  unmatched_symbols_map <- unmatched_symbols_map[is.na(unmatched_symbols_map$Symbol), ]

  # copy the fixed symbols back to the main data structure df2
  df2$response_component <- genes_map$Symbol

  # get rid of genes that had no valid symbol.
  df2 <- df2[!is.na(df2$response_component), ]

  summary_df <- add_to_summary(summary_df,
                               "Number of valid gene symbols" ,
                               length(unique(df2$response_component)))

} else if (response_type == "CELLTYPE_FREQUENCY") {

  # ctf_match will have NA values where no match found
  ctf_match <- match(df2$response_component_original, ctf_fixes$original_annotation)
  # Test if any rows in cell-type mapping sheet are unused.
  w <- !(1:length(ctf_fixes$original_annotation) %in% unique(ctf_match))
  if(any(w)) {
#    print("not all cell-type mappings used:")
#    print(sapply(ctf_fixes$original_annotation[w], function(x) {which(ctf_fixes$original_annotation == x) + 1}))
  } else {
    print("all cell-type mappings used")
  }

  if (any(is.na(ctf_match))) {
    w <- which(is.na(ctf_match))
    print(paste("row", df2$uniq_obs_id[w], ", no mapping for cell type",
                df2$response_component_original[w]))
    missing_mappings <- data.frame(row = df2$uniq_obs_id[w], label = df2$response_component_original[w])
    write.table(missing_mappings,
                file = logfile_path(logdir, base_filename, "missing_mappings.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    write.table(unique(sort(df2$response_component_original[w])),
                file = logfile_path(logdir, base_filename, "missing_mappings_unique.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

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
                file = logfile_path(logdir, base_filename, "no_cell_type_match.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    write.table(sort(failed_matches_uniq),
                file = logfile_path(logdir, base_filename, "no_cell_type_match_sorted.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
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
  df2$response_component_id   <- ctf_fixes$CL_ID[ctf_match]
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

df2 <- cSplit(df2, "exposure_material_id", sep = ";", direction = "long")

if(any(df2$exposure_material_id == "")) {
  # Do not expect this
  stop("missing exposure_material_id")
}
# In case we forget to stop...
df2 <- df2[df2$exposure_material_id != "", ]

write_unique_list(df2$exposure_material_id, logdir, base_filename, "exposure_material_id", do_split = FALSE)

end_cnt <- nrow(df2)
summary_df <- add_to_summary(summary_df, "Split 2, on exposure material_id, rows added", end_cnt - start_cnt)
summary_df <- add_to_summary(summary_df, "Split 2, on exposure material_id, total rows", end_cnt)
summary_df <- add_to_summary(summary_df, "unique ID codes in exposure_material_id", length(unique(df2$exposure_material_id)))
start_cnt <- end_cnt

# The cSplit() calls produce a data.table of data.frame rows.  Coerce back to just data.frame.
df2 <- as.data.frame(df2, stringsAsFactors = FALSE)



###############################################################
##### Data splitting (3): target pathogens (VACCINE ONLY) #####
###############################################################

if (exposure_type == "VACCINE") {

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
  # influenza rows don't hae taxonid yet, supply a dummy so split does not remove rows
  df2 <- cSplit(df2, "target_pathogen_taxonid", sep = ";", direction = "long")

  # class(df2)  # "data.table" "data.frame"
  df2 <- as.data.frame(df2)

  end_cnt <- nrow(df2)
  summary_df <- add_to_summary(summary_df, "Split 3, on pathogens, rows added", end_cnt - start_cnt)
  summary_df <- add_to_summary(summary_df, "Split 3, on pathogens, total rows", end_cnt)
  start_cnt <- end_cnt

  # in Infection template, have allowed "(name)" exposure_material_id column.  Remove.
  codes <- sub(" .*$", "", df2$target_pathogen_taxonid)
  codes <- sub("ncbi_taxid:", "", codes)
  s <- sapply(as.integer(codes), is.integer)
  if (!all(s)) {
    print("target_pathogen_taxonid: unexpected namespace tag")
  }
  df2$target_pathogen_taxonid <- codes

  # Write unique list of target pathogens after substitutions
  text  <- strsplit(df2$target_pathogen, ";") # returns a list
  text <- unique(trimws(unlist(text)))

  # lowercase first letter of each pathogen
#  substr(text, 1, 1) <- tolower(substr(text, 1,1))
  write_unique_list(text, logdir, base_filename, "target_pathogens_after_fixes")
  summary_df <- add_to_summary(summary_df, "target_pathogens_after_fixes", length(text))

  write_unique_list(df2$target_pathogen_taxonid, logdir, base_filename, "target_pathogen_taxonid")
  summary_df <- add_to_summary(summary_df, "target_pathogen_taxonid", length(unique(df2$target_pathogen_taxonid)))
}


###############################################################
##### Data splitting (4): tissue_type_term_id             #####
###############################################################

##### split ######
# influenza rows don't hae taxonid yet, supply a dummy so split does not remove rows
df2 <- cSplit(df2, "tissue_type_term_id", sep = ";", direction = "long")

# class(df2)  # "data.table" "data.frame"
df2 <- as.data.frame(df2)

end_cnt <- nrow(df2)
summary_df <- add_to_summary(summary_df, "Split 4, on tissues, rows added", end_cnt - start_cnt)
summary_df <- add_to_summary(summary_df, "Split 4, on tissues, total rows", end_cnt)
start_cnt <- end_cnt

df2$tissue_type_term_id <- sub(" .*$", "", df2$tissue_type_term_id)

# Write unique list of observed tissues
text  <- strsplit(df2$tissue_type, ";") # returns a list
text <- unique(trimws(unlist(text)))

write_unique_list(text, logdir, base_filename, "tissues_observed")
summary_df <- add_to_summary(summary_df, "tissues_observed", length(text))

write_unique_list(df2$tissue_type_term_id, logdir, base_filename, "tissue_type_term_id")
summary_df <- add_to_summary(summary_df, "tissue_type_term_id", length(unique(df2$tissue_type_term_id)))


####################################################################
#### Get counts of each response component for e.g. word cloud #####
####################################################################

response_df <- unique(df2[ , c("response_component", "uniq_obs_id")])
response_df <- as.data.frame(table(response_df$response_component))
colnames(response_df) <- c("response_component", "count")
response_df <- response_df[order(response_df$count, decreasing = TRUE), ]

# Write out counts in tab-delimited format
# Use .txt extension rather than .tsv to thwart Excel auto-change of symbols
write.table(response_df,
            file = logfile_path(logdir, base_filename, "response_component_counts.txt"),
            sep = "\t", row.names = FALSE)


##########################################################
####  Collect publication titles, dates and abstracts ####
##########################################################

pmids_uniq <- unique(df2$publication_reference_id)

summary_df <- add_to_summary(summary_df, "Unique PMIDs", length(pmids_uniq))
# Get user-entered year for each pmid
# FIXME - need to catch if user did not enter year.
print_pub_year <- df2$publication_date[match(pmids_uniq, df2$publication_reference_id)]
# FIXME - need to check date format valid
print_pub_year <- stri_sub(print_pub_year, 1, 4)

# pmids <- "16571413" # for testing
# pmids <- "24336226" # for testing
# pmids <- "23594957" # for testing
if (RENEW_PMIDS || !file.exists(pmid_file)) {
  td <- mapply(pmid_to_title_easy, pmids_uniq, print_pub_year, SIMPLIFY = FALSE)
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
# FIXME - this is also logged below at comment "original template rows deleted"
which(!(names(signatures_uniq) %in% uniqIDs))
signatures_uniq <- signatures_uniq[names(signatures_uniq) %in% uniqIDs]


# create data structures needed for mSigDB
if (response_type == "GENE" && CREATE_MSIGDB) {
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

  if (response_type == "GENE") {
    base_row$response_component <- paste(unique(df2tmp$response_component), collapse = "; ")
    resp_components_annotated[[i]] <- c(response_rowname, response_description, unique(df2tmp$response_component))

    w <- which(titles_and_dates_df$pmid == sub("pmid:", "", base_row$publication_reference_id))
    if (length(w) != 1) {
      stop(paste("unexpected PMID result: uniqID = ", uniqIDs[i], ", pmid = ", base_row$publication_reference_id))
    }
    titles_and_dates_row <- titles_and_dates_df[w, ]

    if(CREATE_MSIGDB) {
      # Create mSigDB submission row
      msigdb_df <- msigdb_process_row(msigdb_empty, df2tmp)
      msigdb_list[[i]] <- msigdb_df
    }

  } else if (response_type == "CELLTYPE_FREQUENCY") {
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
  if(exposure_type == "VACCINE") {
    base_row$target_pathogen_taxonid   <- paste(unique(df2tmp$target_pathogen_taxonid), collapse = "; ")
  }
  base_row$exposure_material_id <- paste(unique(df2tmp$exposure_material_id), collapse = "; ")
  base_row$tissue_type_term_id <- paste(unique(df2tmp$tissue_type_term_id), collapse = "; ")


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
if (response_type == "GENE" && CREATE_MSIGDB) {
  summary_df <- write_msigdb_submission(msigdb_list, summary_df)
}

# Delete unneeded columns from the recreated template
if (response_type == "GENE") {
  # do nothing
} else if (response_type == "CELLTYPE_FREQUENCY") {
  del_cols <- c("proterm_and_extra", "fully_qualified_response_component")
  recreated_template_df <- recreated_template_df[!colnames(recreated_template_df) %in% del_cols]
}

del_cols <- c("submission_name", "submission_date", "template_name", "short_comment", "process_note")
recreated_template_df <- recreated_template_df[!colnames(recreated_template_df) %in% del_cols]

# Set that first column name back to blank
colnames(recreated_template_df)[1] <- ""
# Write out the recreated upload template in tab-delimited format
write.table(recreated_template_df,
            file = logfile_path(logdir, base_filename, "recreated_template.tsv"),
            sep = "\t", row.names = FALSE)


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
if (exposure_type == "VACCINE") {
   del_cols <- c("exposure_material", "target_pathogen", "short_comment", "process_note")
} else {
  del_cols <- c("exposure_material", "short_comment", "process_note")
}

df2 <- df2[!colnames(df2) %in% del_cols]
header_rows <- header_rows[!colnames(header_rows) %in% del_cols]

if (response_type == "GENE") {
  write_submission_template(df2, header_rows, template_name, titles_and_dates_df,
                            resp_components_collected, unmatched_symbols_map,
                            response_type, exposure_type, project)
} else if (response_type == "CELLTYPE_FREQUENCY") {
  write_submission_template(df2, header_rows, template_name, titles_and_dates_df,
                            resp_components_full_sig, unmatched_symbols_map = NULL,
                            response_type, exposure_type, project)
}

###########################################
##### Write out various summary files #####
###########################################

summary_df <- add_to_summary(summary_df, "total signatures", length(resp_components_collected))

dd0 <- sum(duplicated(resp_components_collected))
summary_df <- add_to_summary(summary_df, "duplicate signatures (response components only)", dd0)
dd0 <- length(unique(resp_components_collected))
summary_df <- add_to_summary(summary_df, "unique signatures (response components only)", dd0)

if (response_type == "GENE") {
  # do nothing
} else if (response_type == "CELLTYPE_FREQUENCY") {
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
rc_cnt_by_pmid <- lapply(pmids_uniq, function(pmid) {
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
rc_cnt_by_pmid_uniq <- resp_comps_by_pmid(resp_components_annotated, pmids_uniq)
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
write_final_summary(all_response_types, exposure_type, logdir)


##########################
##### Generate Plots #####
##########################

# Plot a histogram of unique response components per PMID
# outfile <- paste(base_filename, "UniqueCountByPMID.png", sep = "-")
# png(outfile)
# if (response_type == "GENE") {
#   p <- ggplot(rc_cnt_by_pmid_uniq, aes(count)) + geom_histogram(bins = 10) + scale_x_log10()
#   p +  labs(x = "genes per PMID") + annotation_logticks(sides = "b")
# } else  if (response_type == "CELLTYPE_FREQUENCY") {
#   p <- ggplot(rc_cnt_by_pmid_uniq, aes(count)) + geom_histogram(bins = 10)
#   p +  labs(x = "Cell types per PMID")
# }
# dev.off()

# Plot a histogram of response_components per signature
# png(logfile_path(logdir, base_filename, "unique_count_by_signature.png"))
# if (response_type == "GENE") {
#   p2 <- ggplot(resp_components_cnt_df, aes(count)) + geom_histogram(bins = 10) + scale_x_log10()
#   p2 +  labs(x = "genes per signature") + annotation_logticks(sides = "b") + theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
# } else  if (response_type == "CELLTYPE_FREQUENCY") {
#   p2 <- ggplot(resp_components_cnt_df, aes(count)) + geom_histogram(bins = 9) + scale_x_continuous(breaks = seq(1, 9, 1))
#   p2 +  labs(x = "Cell types per signature") + theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
# }
# dev.off()

#################################
##### Write out sessionInfo #####
#################################

writeLines(capture.output(sessionInfo()),
           logfile_path(logdir, base_filename, "session_info.txt"))

