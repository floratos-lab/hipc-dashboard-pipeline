# FIXME - vaccine_by_year not defined in this file
# If vaccine code is "VO_0000410" ("Pandemrix"), return its single viral component.
# Otherwise, return the components for the specified year, which may also  be coded
getVaccineComponents <- function(year, code) {
  if (code == "VO_0000410") {
    # "pandemrix" only has an entry in the first slot
    vv <- vaccines_by_year[tolower(vaccines_by_year$vaccine.type) == "pandemrix", "name_1"]
  } else {
    w <- grep("name_", colnames(vaccines_by_year))
    vv <- vaccines_by_year[vaccines_by_year$year == year, w]
  }
  vv <- vv[vv != ""]
  return(vv)
}
# getVaccineComponents(2009, "VO_0000410")  # test string

# Return a unique list of pathogens if encounter mixed type of Pandemrix plus Fluvirin.
decodeExposureComponentToPathogen <- function(year, code) {
  # FIXME - fragile, should parse these out
  if (code == "VO_0000046; VO_0000410") {  # Fluvirin; Pandemrix (only occurs for 2009)
    pathogens <- getVaccineComponents("any", "Pandemrix")
    pathogens <- c(pathogens, getVaccineComponents(year, code))  # only uses year
    pathogens <- unique(pathogens)
  }  else {
    pathogens <-  getVaccineComponents(year, code)  # only uses year
  }
  return(pathogens)
}

# Substitute in actual virus components for influenza pathogens only.
# FIXME - if year is not matched in lookup spreadsheet, an undecipherable error occurs.
lookup_vaccine <- function(vaccine_year, exposure_material) {
  # vaccine_year is a comma separated list of one or more years
  # FIXME - note should never actually expect a return value of "unknown" or "",
  #         should catch this as an error
  if (vaccine_year != "") {
    years <- unlist(strsplit(vaccine_year, ","))
    years <- trimws(years)

    if (length(years) == 1) {
      vc <- decodeExposureComponentToPathogen(years, exposure_material)
    } else {
      vc <- sapply(years, decodeExposureComponentToPathogen, exposure_material)
      vc <- unique(unlist(vc))
    }
    if (length(vc) == 0) {
      vc <- "unknown"
    }
    target_pathogen <- paste(vc, collapse = "; ")
  } else if (exposure_material == "VO_0000410") {  # Pandemrix
    vc <- decodeExposureComponentToPathogen("none", exposure_material)
    if (length(vc) == 0) {
      vc <- "unknown"
    }
    # FIXME - Pandemrix only has one component, don't need to collapse...
    target_pathogen <- paste(vc, collapse = "; ")
  } else {
    target_pathogen <- ""
  }
  return(target_pathogen)
}
