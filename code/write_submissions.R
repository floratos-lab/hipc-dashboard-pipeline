generate_observation_summary <- function(sheet_type,
                                         exposure_type,
                                         joining_preposition,
                                         age_string,
                                         tissue_cnt,
                                         exposure_cnt,
                                         pathogen_cnt_os,
                                         at_time_point_phrase = "at <time_point> <time_point_units> from <baseline_time_event>, ",
                                         use_additional = FALSE) {
  # generate a wild-card template, add index if multiple to match new columns in data
  gen_phrase <- function(cnt, phrase) {
    if (cnt == 1) {
      paste0("<", phrase, ">")
    } else {
      new_cols <- paste0("<", phrase, "_", 1:cnt, ">")
      paste0(paste(new_cols[-cnt], collapse = ", "), " and/or ", new_cols[cnt])
    }
  }

  obs_summary <- paste("In", gen_phrase(tissue_cnt, "tissue_type_term_id"))

  if (sheet_type == "GENE") {
    obs_summary <- paste0(obs_summary, ", <response_component> <response_behavior_type>")
  } else if (sheet_type == "CELLTYPE_FREQUENCY") {
    obs_summary <- paste0(obs_summary, ", the frequency of cell type <response_component_id> <proterm_and_extra>")
  }
  obs_summary <- paste0(
    obs_summary, " was <response_behavior> ", at_time_point_phrase,
    joining_preposition, " <comparison>, in cohort <cohort>", age_string,
    ", after exposure to ", gen_phrase(exposure_cnt, "exposure_material_id")
  )
  # Note - pathogen_cnt_os is set to zero when don't want to display pathogens
  if (exposure_type == "VACCINE" &&
    pathogen_cnt_os > 0 && pathogen_cnt_os <= 3) {
    obs_summary <- paste(
      obs_summary, "targeting",
      gen_phrase(pathogen_cnt_os, "target_pathogen_taxonid")
    )
    if (use_additional) {
      obs_summary <- paste(obs_summary, "as well as <additional_exposure_material>")
    }
  }

  return(obs_summary)
}

choose_joining_preposition <- function(response_behavior) {
  if (grepl("^up|^down|enriched", response_behavior, ignore.case = TRUE)) {
    return("for comparison")
  } else if (grepl("correlated|associated", response_behavior, ignore.case = TRUE)) {
    return("with")
  } else if (grepl("predictive", response_behavior, ignore.case = TRUE)) {
    return("of")
  } else if (response_behavior == "differentially expressed") {
    return("for")
  } else {
    warning("choose_joining_preposition(): unexpected response_behavior:", response_behavior)
    return("for")
  }
}

# Write out the submission templates
write_submission_template <- function(df2_cpy, header_rows, release_files, csv_files, template_name, titles_and_dates_df,
                                      resp_components, unmatched_symbols_map = NULL,
                                      sheet_type, exposure_type, project) {
  # create template directory and its file subdirectory
  double_path <- paste(release_files, template_name, "files", sep = "/")
  if (!dir.exists(double_path)) {
    dir.create(double_path, recursive = TRUE)
  }

  cv <- data.frame() # initialize
  uniq_ids_local <- as.character(unique(df2_cpy$sig_row_id))
  print(paste("uniq_ids match:", all(names(resp_components) == uniq_ids_local), collapse = ", "))

  # For each submission, i.e. sig_row_id (one signature)
  for (i in 1:length(uniq_ids_local)) {
    dftmp <- df2_cpy[df2_cpy$sig_row_id == uniq_ids_local[i], ]

    # Despite best efforts, exposure_material_id was being treated as a factor.
    # Remember the data is still split at this point
    uniqExpMat <- unique(as.character(dftmp$exposure_material_id))
    if (exposure_type == "VACCINE") {
      uniqPathogens <- unique(as.character(dftmp$target_pathogen_taxonid))
    }
    uniq_tissue_vals <- unique(as.character(dftmp$tissue_type_term_id))

    # For one submission, all values in these columns are the same, so take first.
    pmid_local <- dftmp$publication_reference_id[1]
    submission_identifier <- dftmp$sig_subm_id[1]
    joining_preposition <- choose_joining_preposition(dftmp$response_behavior[1])

    age_min <- dftmp$age_min[1]
    age_max <- dftmp$age_max[1]
    age_units <- dftmp$age_units[1]

    # get titles and dates for matching PMID
    w <- which(titles_and_dates_df[, "pmid"] == pmid_local)
    if (length(w) != 1) {
      stop("missing PMID result: uniqID = ", uniq_ids_local[i], ", pmid = ", pmid_local)
    }

    dashboard_title <- titles_and_dates_df[w, "dashboard_title"]
    header_rows$publication_reference_url[6] <- dashboard_title

    submission_name <- paste(template_name, pmid_local, submission_identifier, sep = "_")
    dftmp$submission_name <- submission_name
    dftmp$template_name <- submission_name
    dftmp$submission_date <- titles_and_dates_df[w, "date"]

    # Reattach the header rows to each submission in turn
    dftmp <- rbind(header_rows, dftmp)
    colnames(dftmp)[1] <- "" # get rid of "X." or other name in first column

    # Now that header is reattached, add any needed new columns
    # Add a new column to hold the link to the signature file with all response component
    signatureFilename <- paste0(paste(template_name, "sig", pmid_local,
      submission_identifier,
      sep = "_"
    ), ".txt")

    dftmp$signature_file <- signatureFilename

    # customized headers by type
    if (sheet_type == "GENE") {
      dftmp$signature_file[1:6] <- c("", "file", "observed", "", "", "response component (genes) file")

      completeSignatureFilename <- paste0(paste(template_name, "sig_complete", pmid_local,
        submission_identifier,
        sep = "_"
      ), ".txt")
      dftmp$signature_file_complete <- completeSignatureFilename
      dftmp$signature_file_complete[1:6] <- c("", "file", "observed", "", "", " response components including non-HGNC genes")
    } else if (sheet_type == "CELLTYPE_FREQUENCY") {
      dftmp$signature_file[1:6] <- c("", "file", "observed", "", "", "response component (cell types) file")
    }

    # Add new columns for multiple exposure materials
    # FIXME - same code repeated three times in following sections, should make function
    exposure_cnt <- length(uniqExpMat)
    if (exposure_cnt > 1) {
      newCols <- paste("exposure_material_id", 1:exposure_cnt, sep = "_")
      dftmp[newCols] <- ""

      for (j in 1:exposure_cnt) {
        if (exposure_type == "VACCINE") {
          dftmp[1:6, newCols[j]] <- c("vaccine", "", "vaccine", "", "", paste("exposure material", j))
        } else if (exposure_type == "INFECTION") {
          dftmp[1:6, newCols[j]] <- c("pathogen", "", "pathogen", "", "", paste("exposure material", j))
        } else {
          stop("unknown exposure_type")
        }
        dftmp[7:nrow(dftmp), newCols[j]] <- uniqExpMat[j]
      }

      # Remove the original column if now has multiple children
      w <- which(colnames(dftmp) == "exposure_material_id")
      if (length(w) != 1) {
        stop(paste("exposure_material column name mismatch", colnames(dftmp)[w]))
      } else {
        dftmp <- dftmp[, -w]
      }
    }

    if (exposure_type == "VACCINE") {
      # Add new columns for multiple pathogens
      pathogen_cnt <- length(uniqPathogens)
      if (pathogen_cnt > 1) {
        newCols <- paste("target_pathogen_taxonid", 1:pathogen_cnt, sep = "_") # new column names
        dftmp[newCols] <- ""

        for (j in 1:pathogen_cnt) {
          dftmp[1:6, newCols[j]] <- c("pathogen", "", "pathogen", "", "", paste("target pathogen", j))
          dftmp[7:nrow(dftmp), newCols[j]] <- uniqPathogens[j]
        }
        # Remove the original column
        w <- which(colnames(dftmp) == "target_pathogen_taxonid")
        if (length(w) != 1) {
          stop(paste("target_pathogen column name mismatch", colnames(dftmp)[w]))
        }
        dftmp <- dftmp[, -w]
      }
    }

    # Add new columns for multiple tissues
    tissue_cnt <- length(uniq_tissue_vals)
    if (tissue_cnt > 1) {
      newCols <- paste("tissue_type_term_id", 1:tissue_cnt, sep = "_")
      dftmp[newCols] <- ""

      for (j in 1:tissue_cnt) {
        dftmp[1:6, newCols[j]] <- c("cell_subset", "", "tissue", "", "", paste("tissue (Cell Ontology)", j))
        dftmp[7:nrow(dftmp), newCols[j]] <- uniq_tissue_vals[j]
      }

      # Remove the original column if now has multiple children
      w <- which(colnames(dftmp) == "tissue_type_term_id")
      if (length(w) != 1) {
        stop(paste("tissue_type_term_id column name mismatch", colnames(dftmp)[w]))
      } else {
        dftmp <- dftmp[, -w]
      }
    }

    # since all exposure materials and pathogens are now in each row,
    # eliminate complete duplicates
    dftmp <- unique(dftmp)

    # Save the response component list to a separate file for inclusion with submission.
    # Enter the filename into the corresponding template
    # There will be a separate signature file for each uniq_id within a submission (PMID)

    # character, not integer, addressing for lists
    uil <- uniq_ids_local[i]

    # Write out list of submission's response components to file
    # Note that only one response_component per line, so don't worry about separator character.
    # Use suffix ".txt" to prevent automatic Excel conversions of gene names
    write.table(resp_components[[uil]],
      file = paste(release_files, template_name, "files", signatureFilename, sep = "/"),
      row.names = FALSE, col.names = FALSE, quote = FALSE
    )

    if (sheet_type == "GENE") {
      # write the complete signature with all genes, including non-HGNC
      ugs_map <- unmatched_symbols_map[unmatched_symbols_map$sig_row_id == uniq_ids_local[i], "alias"]
      # ugs_map <- sort(ugs_map$alias)
      completeGenes <- c(resp_components[[uil]], ugs_map)
      write.table(completeGenes,
        file = paste(release_files, template_name, "files", completeSignatureFilename, sep = "/"),
        row.names = FALSE, col.names = FALSE, quote = FALSE
      )
    } else if (sheet_type == "CELLTYPE_FREQUENCY") {
      # do nothing
    }

    # Write out the actual submission template in tab-delimited format
    write.table(dftmp,
      file = paste(release_files, template_name, paste0(submission_name, ".txt"), sep = "/"),
      sep = "\t", row.names = FALSE, qmethod = "double"
    )

    if (exposure_type == "VACCINE") {
      # FIXME - vaccine_VO_has_pathogens is a global variable defined in main routine
      if (any(uniqExpMat %in% vaccine_VO_has_pathogens)) {
        pathogen_cnt_os <- 0 # do not display pathogens in observation summary
      } else {
        pathogen_cnt_os <- pathogen_cnt
      }
    } else {
      pathogen_cnt_os <- 0
    }

    # We have to break the usual template wild-card rules and hard code
    # the age values because of the complexity of age combinations.
    age_string <- ""
    if (age_min != "" && age_max != "") {
      if (age_min == age_max) {
        age_string <- age_min
      } else {
        age_string <- paste0(age_min, "-", age_max)
      }
    } else if (age_min != "" && age_max == "") {
      age_string <- paste0(age_min, "+")
    } else if (age_min == "" && age_max != "") {
      age_string <- paste0("<=", age_max)
    }
    if (age_string != "") {
      age_string <- paste(", for age group", age_string, age_units, "old")
      if (dftmp$cohort[1] != "") {
        age_string <- paste0(age_string, ",")
      }
    }

    at_time_point_phrase <-
      "at <time_point> <time_point_units> from <baseline_time_event>, "

    # dtmp is expected to be a data.frame of multiple rows
    # starting row 7 (actual data), many thins are supposed to be the same,
    # at least for time_point and baseline_time_event
    if (trimws(dftmp$baseline_time_event[7]) == "") {
      at_time_point_phrase <- ", "
    } else {
      time_point <- trimws(dftmp$time_point[7])
      if (time_point == "" || time_point == "0") {
        at_time_point_phrase <- "at <baseline_time_event>, "
      }
    }

    observation_summary <- generate_observation_summary(
      sheet_type,
      exposure_type,
      joining_preposition,
      age_string,
      tissue_cnt,
      exposure_cnt,
      pathogen_cnt_os,
      at_time_point_phrase,
      exposure_type == "VACCINE" && dftmp$additional_exposure_material[1] != ""
    )

    if (sheet_type == "GENE") {
      subm_type <- "gene"
    } else if (sheet_type == "CELLTYPE_FREQUENCY") {
      subm_type <- "cell-type freq"
    }

    # generate dashboard-CV-per-template
    cv <- rbind(cv, data.frame(
      observation_tier = 1,
      template_name = submission_name,
      observation_summary = observation_summary,
      story_title = "",
      submission_name = submission_name,
      submission_description = dashboard_title,
      project = project,
      submission_story = "false",
      submission_story_rank = 0,
      submission_center = "HIPC-II Signatures Project",
      principal_investigator = "HIPC-II Sigs: Steven H. Kleinstein, Ph.D.",
      pmid = pmid_local,
      submission_type = subm_type,
      stringsAsFactors = FALSE
    ))
  }

  # Write out the actual dashboard-CV-per-template template in tab-delimited format
  write.table(cv,
    paste(release_files, template_name, paste0(template_name, "-CV-per-template.txt"), sep = "/"),
    sep = "\t", row.names = FALSE
  )
}
