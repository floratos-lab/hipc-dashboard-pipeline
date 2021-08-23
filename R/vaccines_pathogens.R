getVaccineComponents <- function(year) {
  w <- grep("name_", colnames(vaccines_by_year))
  vby <- vaccines_by_year[vaccines_by_year$year == year, w]
  vby <- vby[vby != ""]
  return(vby)
}

# Substitute in actual virus components for influenza pathogens only.
# FIXME - if year is not matched in lookup spreadsheet, an undecipherable error occurs.
lookup_vaccine <- function(vaccine_year) {
  # vaccine_year is a comma separated list of one or more years
  # FIXME - note should never actually expect a return value of "unknown" or "",
  #         should catch this as an error
  if (vaccine_year != "") {
    years <- unlist(strsplit(vaccine_year, ","))
    years <- trimws(years)

    vc <- sapply(years, getVaccineComponents)
    vc <- unique(unlist(vc))

    if (length(vc) == 0) {
      vc <- "unknown"
    }
    target_pathogen <- paste(vc, collapse = "; ")
  }  else {
    target_pathogen <- ""
  }
  return(target_pathogen)
}
