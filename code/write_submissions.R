
submission_center      <- "HIPC-II Signatures Project"
principal_investigator <- "HIPC-II Sigs: Steven H. Kleinstein, Ph.D."

# NOTE - uses "GENE" and "CELLTYPE_FREQUENCY" definitions from main routine
generate_observation_summary <- function(sheet_type,
                                         exposure_type,
                                         joining_preposition,
                                         age_string,
                                         tissue_cnt,
                                         exposure_cnt,
                                         pathogen_cnt,
                                         use_additional = FALSE) {

  # generate a wild-card template, add index if multiple to match new columns in data
  gen_phrase <- function(cnt, phrase) {
    if(cnt == 1) {
      os <- paste0("<", phrase, ">")
    } else {
      newCols <- sapply(1:cnt, function(cntr) {
        paste0("<", phrase, "_", cntr, ">")
      })
      os <- paste(newCols[1:(cnt-1)], collapse = ", ")
      os <- paste(os, newCols[cnt], sep = " and/or ")
    }
    return(os)
  }

  obs_summary <- "In"
  obs_summary <- paste(obs_summary, gen_phrase(tissue_cnt, "tissue_type_term_id"))

  if (sheet_type == "GENE") {
     obs_summary     <- paste0(obs_summary, ", <response_component> <response_behavior_type>")
  } else if (sheet_type == "CELLTYPE_FREQUENCY") {
    # FIXME - the responses may not be just frequency.  There can also be activation state.
    #         This could be implemented using <response_behavior_type> to specify actual.
    obs_summary     <- paste0(obs_summary,", <response_component_id> <proterm_and_extra> frequency")
  }
  obs_summary   <- paste(obs_summary, "was <response_behavior>")
  obs_summary   <- paste(obs_summary, "at <time_point> <time_point_units> from <baseline_time_event>")
  obs_summary   <- paste(obs_summary, joining_preposition)
  obs_summary   <- paste(obs_summary, "<comparison> in")
  obs_summary   <- paste(obs_summary, "cohort", "<cohort>", age_string)
  obs_summary <- paste(obs_summary, "after exposure to")
  obs_summary <- paste(obs_summary, gen_phrase(exposure_cnt, "exposure_material_id"))
  # Note - pathogen_cnt is set to zero when don't want to display pathogens
  if(exposure_type == "VACCINE") {
    if (pathogen_cnt > 0 && pathogen_cnt <= 3) {
      obs_summary <- paste(obs_summary, "targeting")
      obs_summary <- paste(obs_summary, gen_phrase(pathogen_cnt, "target_pathogen_taxonid"))
      if(use_additional) {
        obs_summary <- paste(obs_summary, "as well as <additional_exposure_material>")
      }
    }
  }

  return(obs_summary)
}

choose_joining_preposition <- function(response_behavior) {
  if(grepl("^up|^down|enriched", response_behavior, ignore.case = TRUE)) {
    prep <- "for comparison"
  } else if (grepl("correlated|associated", response_behavior, ignore.case = TRUE)) {
    prep <- "with"
  } else if (grepl("predictive",  response_behavior, ignore.case = TRUE)) {
    prep <- "of"
  } else {
    print(paste("choose_joining_preposition(): unexpected response_behavior:", response_behavior))
    prep <- "for"
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

  # A place to keep csv format copies of the templates for easy inspection in Excel
  template_name_csv <- paste(template_name, "csv", sep = "_")
  double_path <- paste(csv_files, template_name_csv, "files", sep = "/")
  if (!dir.exists(double_path)) {
    dir.create(double_path, recursive = TRUE)
  }

  cv <- data.frame()  # initialize
  uniq_ids_local <- as.character(unique(df2_cpy$sig_row_id))
  print(paste("uniq_ids match:", all(names(resp_components) == uniq_ids_local), collapse = ", "))

  # For each submission, i.e. sig_row_id (one signature)
  for (i in 1:length(uniq_ids_local)) {
    dftmp <- df2_cpy[df2_cpy$sig_row_id == uniq_ids_local[i], ]

    # Despite best efforts, exposure_material_id was being treated as a factor.
    # Remember the data is still split at this point
    uniqExpMat <- unique(as.character(dftmp$exposure_material_id))
    if(exposure_type == "VACCINE") {
      uniqPathogens <- unique(as.character(dftmp$target_pathogen_taxonid))
    }
    uniq_tissue_vals <- unique(as.character(dftmp$tissue_type_term_id))

    # For one submission, all values in these columns are the same, so take first.
    pmid_local <- dftmp$publication_reference_id[1]
    submission_identifier <- dftmp$sig_subm_id[1]
    joining_preposition <- choose_joining_preposition(dftmp$response_behavior[1])

    age_min   <- dftmp$age_min[1]
    age_max   <- dftmp$age_max[1]
    age_units <- dftmp$age_units[1]

    if(exposure_type == "VACCINE") {
      use_additional <- dftmp$additional_exposure_material[1] != ""
    } else {
      use_additional <- FALSE
    }

    # get titles and dates for matching PMID
    w <- which(titles_and_dates_df[, "pmid"] == pmid_local)
    if (length(w) != 1) {
      stop("missing PMID result: uniqID = ", uniq_ids_local[i], ", pmid = ", pmid_local)
    }

    dashboard_title <- titles_and_dates_df[w, "dashboard_title"]
    header_rows$publication_reference_url[6] <- dashboard_title

    submission_name <- paste(template_name, pmid_local, submission_identifier, sep = "_")
    dftmp$submission_name <- submission_name
    dftmp$template_name   <- submission_name
    dftmp$submission_date <- titles_and_dates_df[w, "date"]

    # Reattach the header rows to each submission in turn
    dftmp <- rbind(header_rows, dftmp)
    colnames(dftmp)[1] <- ""  # get rid of "X." or other name in first column

    # Now that header is reattached, add any needed new columns
    # Add a new column to hold the link to the signature file with all response component
    signatureFilename <- paste0(paste(template_name, "sig", pmid_local,
                                      submission_identifier, sep = "_"), ".txt")

    dftmp$signature_file <- signatureFilename

    # customized headers by type
    if (sheet_type == "GENE") {
      dftmp$signature_file[1:6] <- c("", "file", "observed", "", "", "response component (genes) file")

      completeSignatureFilename <- paste0(paste(template_name, "sig_complete", pmid_local,
                                                submission_identifier, sep = "_"), ".txt")
      dftmp$signature_file_complete <- completeSignatureFilename
      dftmp$signature_file_complete[1:6] <- c("", "file", "observed", "", "", " response components including non-HGNC genes")
    } else if (sheet_type == "CELLTYPE_FREQUENCY") {
      dftmp$signature_file[1:6] <- c("", "file", "observed", "", "", "response component (cell types) file")
    }

    # Add new columns for multiple exposure materials
    expMatCnt <- length(uniqExpMat)
    if(expMatCnt > 1) {
      newCols <- paste("exposure_material_id", 1:expMatCnt, sep = "_")
      dftmp[newCols] <- ""

      for (j in 1:expMatCnt) {
        if(exposure_type == "VACCINE") {
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
      if(length(w) != 1) {
        stop(paste("exposure_material column name mismatch", colnames(dftmp)[w]))
      } else {
        dftmp <- dftmp[, -w]
      }
    }

    if(exposure_type == "VACCINE") {
      # Add new columns for multiple pathogens
      pathogen_cnt <- length(uniqPathogens)
      if(pathogen_cnt > 1) {
        newCols <- paste("target_pathogen_taxonid", 1:pathogen_cnt, sep = "_")  # new column names
        dftmp[newCols] <- ""

        for (j in 1:pathogen_cnt) {
          dftmp[1:6, newCols[j]] <- c("pathogen", "", "pathogen", "", "", paste("target pathogen", j))
          dftmp[7:nrow(dftmp), newCols[j]] <- uniqPathogens[j]
        }
        # Remove the original column
        w <- which(colnames(dftmp) == "target_pathogen_taxonid")
        if(length(w) != 1) {
          stop(paste("target_pathogen column name mismatch", colnames(dftmp)[w]))
        }
        dftmp <- dftmp[, -w]
      }
    }

    # Add new columns for multiple tissues
    tissue_cnt <- length(uniq_tissue_vals)
    if(tissue_cnt > 1) {
      newCols <- paste("tissue_type_term_id", 1:tissue_cnt, sep = "_")
      dftmp[newCols] <- ""

      for (j in 1:expMatCnt) {
          dftmp[1:6, newCols[j]] <- c("cell subset", "", "tissue", "", "", paste("tissue (Cell Ontology)", j))
        dftmp[7:nrow(dftmp), newCols[j]] <- uniq_tissue_vals[j]
      }

      # Remove the original column if now has multiple children
      w <- which(colnames(dftmp) == "tissue_type_term_id")
      if(length(w) != 1) {
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
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    if (sheet_type == "GENE") {
      # write the complete signature with all genes, including non-HGNC
      ugs_map <- unmatched_symbols_map[unmatched_symbols_map$sig_row_id == uniq_ids_local[i], "alias"]
      # ugs_map <- sort(ugs_map$alias)
      completeGenes <- c(resp_components[[uil]], ugs_map)
      write.table(completeGenes,
                  file = paste(release_files, template_name, "files", completeSignatureFilename, sep = "/"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    } else if (sheet_type == "CELLTYPE_FREQUENCY") {
      # do nothing
    }

    # Write a copy of the response components to the submission CSV directory
    # Only difference is the files end in .csv so can automatically open in Excel
    # For internal use only, not for users as Excel will alter some gene symbols.
    write.table(resp_components[[uil]],
                file = paste(csv_files, template_name_csv, "files", signatureFilename, sep = "/"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    if (sheet_type == "GENE") {
      write.table(completeGenes,
                  file = paste(csv_files, template_name_csv, "files", completeSignatureFilename, sep = "/"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    } else if (sheet_type == "CELLTYPE_FREQUENCY") {
      # do nothing
    }

    # Write out the actual submission template in tab-delimited format
    write.table(dftmp,
                file = paste(release_files, template_name, paste0(submission_name, ".txt"), sep = "/"),
                sep = "\t", row.names = FALSE)

    # For convenience, write out CSV version
    write.csv(dftmp,
              file = paste(csv_files, template_name_csv, paste0(submission_name, ".csv"), sep = "/"),
              row.names = FALSE)

    if(exposure_type == "VACCINE") {
      # FIXME - vaccine_VO_has_pathogens is a global variable defined in main routine
      if (any(uniqExpMat %in% vaccine_VO_has_pathogens)) {
        pathogen_cnt_os <- 0  # do not display pathogens in observation summary
      } else {
        pathogen_cnt_os <- pathogen_cnt
      }
    } else {
      pathogen_cnt_os <- 0
    }

    # We have to break the usual template wild-card rules and hard code
    # the age values because of the complexity of age combinations.
    age_string <- ""
    if(age_min != "" && age_max != "") {
      if(age_min == age_max) {
        age_string <- age_min
      } else {
        age_string <- paste0(age_min, "-", age_max)
      }
    } else if(age_min != "" && age_max == "") {
      age_string <- paste0(age_min, "+")
    } else if(age_min == "" && age_max != "") {
      age_string <- paste0("<=", age_max)
    }
    if(age_string != "") {
      if (age_units != "") {
        age_units_rewritten <- ifelse(age_units == "years", "yo",
               ifelse(age_units == "months", "mo",
               ifelse(age_units == "weeks", "wo", age_units)))
        age_string <- paste(age_string, age_units_rewritten)
      }
      if(dftmp$cohort[1] != "") {
        age_string <- paste0(age_string, ",")
      }
    }

    observation_summary <- generate_observation_summary(sheet_type,
                                                        exposure_type,
                                                        joining_preposition,
                                                        age_string,
                                                        tissue_cnt,
                                                        expMatCnt,
                                                        pathogen_cnt_os,
                                                        use_additional)

    if (sheet_type == "GENE") {
      subm_type <- "gene"
    } else if (sheet_type == "CELLTYPE_FREQUENCY") {
      subm_type <- "cell-type freq"
    }

    # generate dashboard-CV-per-template
    cv <- rbind(cv, data.frame(observation_tier = 1,
                               template_name = submission_name,
                               observation_summary = observation_summary,
                               story_title = "",
                               submission_name = submission_name,
                               submission_description = dashboard_title,
                               project = project,
                               submission_story = "false",
                               submission_story_rank = 0,
                               submission_center = submission_center,
                               principal_investigator = principal_investigator,
                               pmid = pmid_local,
                               submission_type = subm_type,
                               stringsAsFactors = FALSE))
  }

  # Write out the actual dashboard-CV-per-template template in tab-delimited format
  write.table(cv,
              paste(release_files, template_name, paste0(template_name, "-CV-per-template.txt"), sep = "/"),
              sep = "\t", row.names = FALSE)

  # Write out the actual dashboard-CV-per-template template in csv format
  write.csv(cv,
            paste(csv_files, template_name_csv, paste0(template_name, "-CV-per-template.csv"), sep = "/"),
            row.names = FALSE)
}