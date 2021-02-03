library(pracma)

# Not really just logging, can be anything needing a parent directory
# and a filename
logfile_path <- function(logdir, base_filename = NULL, filename) {
  if(is.null(base_filename)) {
    outfile <- paste(logdir, filename, sep = "/")
  } else {
    outfile <- paste0(logdir, "/", base_filename, "-", filename)
  }
  return(outfile)
}

add_to_summary <- function(x, text, value) {
  return(rbind(x, data.frame(text = text, value = value)))
}

write_unique_list <- function(inputList, logdir, base_filename, file_description, do_split = FALSE) {
  # Using do_split lets this function work on the original, unsplit data.
  if (do_split) {
    # as.character added because there is a non-printing or otherwise invalid character somewhere in the data...
    val <- strsplit(as.character(inputList),";") # the input column will split to a list of lists
  } else {
    val <- inputList
  }
  # unspool list structure
  write.table(unique(sort(trimws(unlist(val)))),
              file = logfile_path(logdir, base_filename, paste0(file_description, ".txt")),
              sep = "\t", row.names = FALSE, col.names = FALSE)
}

# Format URLs and handle case of multiple URLs
# Due to length restriction of 255 characters on evidence columns,
# can in practice only return at most two concatenated URLs.
# NOTE - There are special cases in this code
format_hipc_urls <- function(url_vector) {
  w <- url_vector != ""
  if(any(w)) {
    s <- strsplit(url_vector[w], "; ")
    # Check if more than 2 URLs are present for any row
    too_long <- sapply(s, function(x) length(x) > 2)
    if(any(too_long)) {
      # NOTE - too_long index not adjusted for header rows
      print(paste("WARNING - more than 2 URLs found at index values:", which(too_long)))
    }
    # Because after strsplit() there can be multiple items per entry in "s",
    # cannot just check in below function if input is null, as will get
    # a length mismatch for > 1 entry.

    # NOTE - the "split" target must be hand-specified for each new case...
    s <- sapply(s, function(x) {
      z <- strsplit(x, split="PMC\\d+|pubmed|pii|pnas\\.\\d+/-|QAD")
      z <- lapply(z, "[", 2)  # get rid of uninteresting part of URL
      z <- sub("/bin/", "", z)
      z <- gsub("^/|/$", "", z) # remove leading and trailing slashes
      z <- sub("NIHMS83712-supplement-", "", z)  # special case, filename has two occurrences of "supplement".
      z <- sub("NIHMS\\d+-", "", z)
      # Form the URLs
      z <- paste0("<a href = '", x, "' target='_blank'>", z, "</a>")
      # concatenate multiple URLs if present
      z <- paste(z, collapse = "; ")
      if(nchar(z) > 255) {
        stop(paste("too many characters in URL", z))
      }
      return(z)
    })
    print(paste("maximum URL string length", max(nchar(s)), sep = ":"))
    url_vector[w] <- s
  }
  return(url_vector)
}

# Unique annotated response component count by PMID
resp_comps_by_pmid <- function(resp_components_annotated, pmids) {

  l3 <- lapply(pmids, function(pmid) {
    # for each PMID, get all the response_components in the list
    l <- lapply(resp_components_annotated, function(ra) {
      if (grepl(pmid, ra[1]) == TRUE) {
        return(ra[-c(1,2)])
      } else {
        return(NULL)
      }
    })
    l <- unique(unlist(l))
    return(list(pmid = pmid, count = length(l)))
  })
  return(as.data.frame(rbindlist(l3)))
}

get_response_type <- function(response_behavior) {
  if(grepl("up|down", response_behavior)) {
    response_type <- UPDOWN
  } else if (grepl("positive|negative|significant", response_behavior)) {
    response_type <- CORRELATION
  } else if (response_behavior == "predictive") {
    response_type <- PREDICTIVESET
  } else if (response_behavior == "enriched") {
    response_type <- ENRICHED
  } else {
    print(paste("unexpected behavior value: ", response_behavior))
    response_type <- ""
  }
  return(response_type)
}

abbreviate_route <- function(route) {
  route <- trimws(route)
  route <- ifelse(route == "", "",
  ifelse(strcmpi(route,"Intramuscular, intradermal, transcutaneous"), "IM, ID, TRANSQ",
  ifelse(strcmpi(route, "Intramuscular injection"), "IM",
  ifelse(strcmpi(route, "i.m."), "IM",
  ifelse(strcmpi(route, "subcutaneous"), "SUBQ",
  ifelse(strcmpi(route, "i.n."), "IN",
  ifelse(strcmpi(route, "PO (oral)"), "PO",
  ifelse(strcmpi(route, "ID (intradermal)"), "ID",
                stop(paste("unknown route encountered:", route))
         ))))))))
  return(route)
}

abbreviate_time_unit <- function(time_unit) {
  time_unit <- ifelse(grepl("hour",  time_unit, ignore.case = TRUE), "H",
               ifelse(grepl("day",   time_unit, ignore.case = TRUE), "D",
               ifelse(grepl("week",  time_unit, ignore.case = TRUE), "W",
               ifelse(grepl("month", time_unit, ignore.case = TRUE), "M",
               ifelse(grepl("year",  time_unit, ignore.case = TRUE), "Y", "")))))
  return(time_unit)
}


