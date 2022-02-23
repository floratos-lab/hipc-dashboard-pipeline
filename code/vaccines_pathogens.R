# FIXME - vaccines_by_year is a global variable

get_vaccine_components <- function(year, type = "taxid") {
  if (type == "taxid") {
    w <- grep("ncbitax_", colnames(vaccines_by_year))
  } else if (type == "name") {
    w <- grep("name_", colnames(vaccines_by_year))
  } else {
    return("")  # FIXME - should be an error
  }
  vby <- vaccines_by_year[vaccines_by_year$year == year, w]
  vby <- vby[vby != ""]
  vby <- vby[!is.na(vby)]
  return(vby)
}

# Substitute in actual virus components for influenza pathogens only.
# FIXME - if year is not matched in lookup spreadsheet, an undecipherable error occurs.
lookup_vaccine <- function(vaccine_year, type) {
  # vaccine_year is a comma separated list of one or more years
  # FIXME - note should never actually expect a return value of "unknown" or "",
  #         should catch this as an error
  if (vaccine_year != "") {
    years <- unlist(strsplit(vaccine_year, ","))
    years <- trimws(years)

    vc <- sapply(years, get_vaccine_components, type)
    vc <- unlist(vc) # a list is returned above if the results are of different length
    vc <- unique(as.vector(vc)) # a matrix is returned above if all results are of same length
    if (length(vc) == 0) {
      vc <- "unknown"
    }
    target_pathogen <- paste(vc, collapse = "; ")
  }  else {
    target_pathogen <- ""
  }
  return(target_pathogen)
}
