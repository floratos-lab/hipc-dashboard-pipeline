# depends on get_response_type() in write_submissions.R

msigdb_intialize <- function() {
  # Columns needed for mSigDB dataframe
  msigdb_cols <- c("standard_name",
                   "organism",
                   "external_details_url",
                   "pmid",
                   "geoid",
                   "exact_source",
                   "chip",
                   "contributor",
                   "contributor_org",
                   "contributor_email",
                   "description_brief",
                   "description_full",
                   "members",
                   "members_updated",
                   "row_key",
                   "short_comment",
                   "process_note")

  msigdb_empty <- data.frame(matrix(data = rep("", length(msigdb_cols)),
                                    nrow=1, ncol=length(msigdb_cols),
                                    dimnames=list(NULL, msigdb_cols)),
                             stringsAsFactors = FALSE)

  return(msigdb_empty)
}


# for MSigDB gene submission
msigdb_description_brief <- function(sheet_type,
                                     tissue_type,
                                     exposure_material,
                                     cohort,
                                     route,
                                     time_point,
                                     time_unit,
                                     response_behavior,
                                     comparison,
                                     comment) {

  if (sheet_type != "GENE") {
    print("brief_description() only supports gene submissions")
    return("")
  }

  obs_summary <- paste("Genes", response_behavior, comparison)
  obs_summary <- paste(obs_summary, "in", tissue_type)
  obs_summary   <- paste(obs_summary, "in", cohort)
  obs_summary   <- paste(obs_summary, "after exposure to", exposure_material)


  time_unit <- abbreviate_time_unit(time_unit)
  obs_summary <- paste(obs_summary, ", time point", paste0(time_point, time_unit))

  if(route != "") {
    obs_summary   <- paste(obs_summary, ", administered", route)
  }

  if(comment != "") {
    obs_summary <- paste(obs_summary, ". Comment: ", comment)
  }

  return(obs_summary)
}

msigdb_standard_name <- function(author,
                                 tissue_type,
                                 exposure_material,
                                 cohort,
                                 route,
                                 time_point,
                                 time_unit,
                                 response_behavior,
                                 comparison,
                                 short_comment) {

  tissue_type <- sub("peripheral blood mononuclear cell", "PBMC", tissue_type, ignore.case = TRUE)

  exposure_material <- trimws(sub("\\.$", "", exposure_material))
  exposure_material <- trimws(gsub("; ", "/", exposure_material))

  exposure_material <- gsub("\\(.*\\)", "", exposure_material)

  exposure_material <- sub("seasonal trivalent inactivated influenza vaccine", "TIV",
                                exposure_material, ignore.case = TRUE)
  exposure_material <- sub("trivalent inactivated vaccine", "TIV",
                                exposure_material, ignore.case = TRUE)
  exposure_material <- sub("inactivated monovalent influenza A", "INACT_MONOV_INFLUENZA_A",
                                exposure_material, ignore.case = TRUE)
  exposure_material <- sub("split-virus vaccine", "",
                                exposure_material, ignore.case = TRUE)
  exposure_material <- trimws(exposure_material)

  response_behavior <- sub("(not previously immunized)", "", response_behavior)

  route <- abbreviate_route(route)
  time_unit <- abbreviate_time_unit(time_unit)

  p <- paste(author, tissue_type, exposure_material, comparison, cohort, route,
             paste0(time_point, time_unit),
             short_comment, response_behavior, sep = "_" )
  p <- gsub("<=", "LTE", p)
  p <- gsub(">=", "GTE", p)
  p <- gsub("<", "LT", p)
  p <- gsub(">", "GT", p)
  p <- gsub("[,.]+", "", p)
  p <- gsub("[ /-]+", "_", p)
  p <- gsub("[\\+]", "PLUS", p)
  p <- gsub("[\\(\\)]+", "", p)
  p <- gsub("_{2,}", "_", p)
  p <- sub("_$", "", p)  # for blank disambiguators
  return(toupper(p))
}

msigdb_process_row <- function(msigdb_empty, df2tmp) {
  msigdb_df <- msigdb_empty
  base_row <- df2tmp[1, ] # get first row for this uniqID

  msigdb_df$standard_name <- msigdb_standard_name(titles_and_dates_row$author,
                                                  base_row$tissue_type,
                                                  base_row$exposure_material,
                                                  base_row$cohort,
                                                  base_row$route,
                                                  base_row$time_point,
                                                  base_row$time_point_unit,
                                                  base_row$response_behavior,
                                                  base_row$comparison,
                                                  base_row$short_comment)

  msigdb_df$organism             <- "Human"
  msigdb_df$external_details_url <- base_row$signature_source_url
  msigdb_df$pmid                 <- base_row$publication_reference
  msigdb_df$geoid                <- ""
  msigdb_df$exact_source         <- base_row$signature_source
  msigdb_df$chip                 <- "Human_Gene_Symbol"
  msigdb_df$contributor          <- "HIPC SIGNATURES"
  msigdb_df$contributor_org      <- "NIAID/HIPC SIGNATURES"
  msigdb_df$contributor_email    <- "feedback@hipc-dashboard.org"

  msigdb_df$description_brief    <- msigdb_description_brief(sheet_type,
                                                             base_row$tissue_type,
                                                             gsub("; *", "/", base_row$exposure_material),
                                                             base_row$cohort,
                                                             base_row$route,
                                                             base_row$time_point,
                                                             base_row$time_point_unit,
                                                             base_row$response_behavior,
                                                             base_row$comparison,
                                                             base_row$comment)

  msigdb_df$description_full <- titles_and_dates_row$abstract
  # for this one column, each row has a separate response component
  msigdb_df$members <- paste(unique(df2tmp$response_component_original), collapse = ",")
  msigdb_df$members_updated <- paste(unique(df2tmp$response_component), collapse = ",")
  msigdb_df$row_key <- base_row$row_key
  msigdb_df$short_comment <- base_row$short_comment
  msigdb_df$process_note <- base_row$process_note
  return(msigdb_df)
}

write_msigdb_submission <- function (msigdb_list, summary_df) {

  msigdb_df <- as.data.frame(rbindlist(msigdb_list))

  # Filter out rows with fewer than 5 genes from mSigDB submission
  nrow(msigdb_df)
  w <- sapply(strsplit(msigdb_df$members, ","), function(x) {length(x) > 4})
  msigdb_df <- msigdb_df[w, ]
  nrow(msigdb_df)
  msigdb_df <- msigdb_df[order(msigdb_df$pmid), ]

  msigdb_standard_name_max_nchar <- max(sapply(msigdb_df$standard_name, nchar))
  summary_df <- add_to_summary(summary_df, "msigdb standard name max nchar", msigdb_standard_name_max_nchar)
  print(paste("msigdb standard name max nchar", msigdb_standard_name_max_nchar))

  # collect all of the rows matching the duplicate entries so can compare
  # FIXME - should explicitly handle case where there are no duplicates, though it works anyway.
  msigdb_dups <- unique(msigdb_df$standard_name[duplicated(msigdb_df$standard_name)])
  s <- sapply(msigdb_dups, function(x) {which(msigdb_df$standard_name == x)})
  ll <- lapply(s, function(x) {msigdb_df[x, ]})
  msigdb_dups <- as.data.frame(rbindlist(ll))

  # Write out the special mSigDB version in Excel format
  write.xlsx(msigdb_df,
             file = logfile_path(logdir, base_filename, "mSigDB_submission.xlsx"),
             sheetName = sheet_name, row.names = FALSE)

  # Write out the duplicated rows from mSigDB version in Excel format
  filep <- logfile_path(logdir, base_filename, "mSigDB_dups.xlsx")
  if(nrow(msigdb_dups) > 0) {
    write.xlsx(msigdb_dups,
               file = filep,
               sheetName = sheet_name, row.names = FALSE)
  } else {
    if(file.exists(filep)) {
      file.remove(filep)
    }
  }
  return(summary_df)
}


