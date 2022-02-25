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

# This was used for mSigDB.  Do not use for regular Dashboard submissions.
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

# check each column of a data frame for the test character
strict_char_check <- function(in_df, testchar) {
  l_out <- lapply(colnames(in_df), function(x) {
    w <- which(sapply(in_df[x], function(y) {grepl(testchar, y)}))
    if(length(w) > 0){
      return(paste(x, paste(w, collapse = ", ")))
    } else {
      return(NULL)
    }
  })
  keep <- !sapply(l_out, is.null)
  if(any(keep)) {
    return(l_out[keep])
  } else {
    return(NULL)
  }
}

# check for lists of response components within one PMID that have a large overlap with one another.
# This is intended to catch the case were one set was accidentally appended to another.
check_response_components_overlap <- function(df2, pmids, min_intersection, min_overlap_fraction, require_different_behaviors, max_hits) {

  final_list <- vector(mode = "list", max_hits)
  cnt <- 1
  for(pmid in pmids) {
    rows_for_pmid <- unique(df2$sig_row_id[df2$publication_reference_id == pmid])
    for(i in 1:length(rows_for_pmid)) {
      rco_i <- df2$response_component_original[df2$sig_row_id == rows_for_pmid[i]]
      rb_i <- df2$response_behavior[df2$sig_row_id == rows_for_pmid[i]][1]
      for(j in min(i + 1, length(rows_for_pmid)):length(rows_for_pmid)) {
        if(i == j) {
          next
        }
        rco_j <- df2$response_component_original[df2$sig_row_id == rows_for_pmid[j]]
        rb_j <- df2$response_behavior[df2$sig_row_id == rows_for_pmid[j]][1]
        inter_set <- intersect(rco_i, rco_j)
        overlap_longer_fraction <- length(inter_set)/max(length(rco_i), length(rco_j))
        overlap_shorter_fraction <- length(inter_set)/min(length(rco_i), length(rco_j))
        if(length(inter_set) < min_intersection) {
          next
        }
        if(overlap_shorter_fraction < min_overlap_fraction) {
          next
        }
        if(require_different_behaviors && (rb_i == rb_j)) {
          next
        }
        final_list[[cnt]] <- data.frame(pmid = pmid,
                                        row1 = rows_for_pmid[i],
                                        row2 = rows_for_pmid[j],
                                        overlap_length = length(inter_set),
                                        set1_count = length(rco_i),
                                        set2_count = length(rco_j),
                                        overlap_longer_fraction = overlap_longer_fraction,
                                        overlap_shorter_fraction = overlap_shorter_fraction,
                                        set1_resp_behav = rb_i,
                                        set2_resp_behav = rb_j
        )
        cnt <- cnt + 1
        if(cnt > max_hits) {
          return(NULL)
        }
      }
    }
  }
  ft <- rbindlist(final_list)
}

