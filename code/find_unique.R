# Returns FALSE if not all types present, TRUE otherwise.
# WARNING - The format of each result file to be read must be specified.
write_joint_summary <- function(response_types, exposure_type, logdir) {
  # check if one file of each response_type is present
  # FIXME - this will not catch rare case of when new file type added in code,
  #         but not yet present on disk.  Should check for each file.
  for(i in 1:length(response_types)) {
    test_filename <- logfile_path(logdir, response_types[i], "exposure_material_id.txt")
    if(!file.exists(test_filename)) {
      print(paste("write_final_summary: file not found", test_filename))
      return(FALSE)
    }
  }

  summarize_response <- function(response_types, type, type_file, header, get_column, logdir) {

    joint_items <- sapply(response_types, function(x) {
      s <- read.table(file = logfile_path(logdir, x, type_file),
                      header = header, stringsAsFactors = FALSE)
      return(s[ , get_column])})

    summary_row <- data.frame(type = type,
                              matrix(data = sapply(joint_items, length), nrow = 1),
                              joint = length(unique(unlist(joint_items))))
    return(summary_row)
  }

  joint_summary <- summarize_response(response_types,
                            type = ifelse(exposure_type == "VACCINE", "vaccines", "pathogen"),
                            type_file = "exposure_material_id.txt",
                            header = FALSE, get_column = "V1",
                            logdir = logdir)

  if(exposure_type == "VACCINE") {
    joint_summary <- rbind(joint_summary,
                           summarize_response(response_types,
                              type = "target_pathogens",
                              type_file = "target_pathogens_after_fixes.txt",
                              header = FALSE, get_column = "V1",
                              logdir = logdir))
  }

  joint_summary <- rbind(joint_summary,
                         summarize_response(response_types,
                             type = "PMIDs",
                             type_file = "response_component_unique_count_by_PMID.csv",
                             header = TRUE, get_column = "pmid",
                             logdir = logdir))

  joint_summary <- rbind(joint_summary,
                         summarize_response(response_types,
                             type = "tissues",
                             type_file = "tissue_type_term_id.txt",
                             header = FALSE, get_column = "V1",
                             logdir = logdir))

  names(joint_summary)[2:(length(response_types)+1)] <- response_types
  rownames(joint_summary) <- joint_summary$type
  joint_summary <- joint_summary[, -1]
  joint_summary <- t(joint_summary)
  write.csv(joint_summary, file = logfile_path(logdir, NULL,
                                  paste(exposure_type, "joint_signatures_count_summary.csv", sep = "_")))
  return(TRUE)
}