# It would be nice to use a direct call like this, then parse the xml return.
# url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=17651872"
# text <- read_html(url)
# But instead, try using easyPubMed.
library(xml2)
library(easyPubMed)

### DATA ###
# Corrections for PMIDs where Pubmed does not return complete publication dates
pmidMonths <- matrix(c("23844129", "Jul",
                       "23658707", "May",
                       "23594957", "Jul",
                       "27529750", "Aug",
                       "25706537", "Feb",
                       "29868000", "May",
                       "30873150", "Feb",
                       "29535712", "Feb",
                       "26148331", "Jul",
                       "28099485", "Jan"),
                     ncol = 2, byrow = TRUE)
pmidMonths <- as.data.frame(pmidMonths)
colnames(pmidMonths) <- c("pmid", "month")

pmidYears <- matrix(c("23594957", "2013"), ncol = 2, byrow = TRUE)
pmidYears <- as.data.frame(pmidYears)
colnames(pmidYears) <- c("pmid", "year")

getPMIDMonth <- function(pmidi) {
  if(pmidi %in% pmidMonths$pmid) {
    rv <- as.character(subset(pmidMonths, pmid == pmidi, select = month))
  } else {
    rv <- ""
  }
  return(rv)
}
# test
# getPMIDMonth("29868000")
# getPMIDMonth("2986800")
# getPMIDMonth("23594957")


getPMIDYear <- function(pmidi) {
  if(pmidi %in% pmidYears$pmid) {
    rv <- as.character(subset(pmidYears, pmid == pmidi, select = year))
  } else {
    rv <- ""
  }
  return(rv)
}
# test
# getPMIDYear("23594957")
# getPMIDYear("2359495")



# test:
# pmid_to_title_easy("16571413")
# pmid_to_title_easy("23594957")

# Returns title and date in two objects.
# Note - There are many ways to get dates from PubMed.
# This way seems to give the actual publication date
# as supplied by the publisher.  However, the
# data is in some cases incomplete, necessitating
# providing some data manually.
# May want to look again at different ways to
# get these dates.

# MEDLINE/PubMed Data Element (Field) Descriptions
# https://www.nlm.nih.gov/bsd/mms/medlineelements.html
# Date of Publication (DP)
# Date of Publication contains the full date on which the issue of the journal was published.
# The standardized format consists of elements for a 4-digit year, a 3-character abbreviated month, and a 1 or 2-digit day.
# Every record does not contain all of these elements;
# the data are taken as they are published in the journal issue,
# with minor alterations by NLM such as abbreviations.
#
# Example:
# DP - 2001 Apr 15
# DP - 2001 Apr
# DP - 2000 Spring
# DP - 2000 Nov-Dec
# DP - 2001


# date types "pubmed", "entrez" and "medline" return properly formated
# full dates.  "medline" is often far later than the real publication date.
# "entrez" and "pubmed" are almost exactly the same as the epub date.
pmid_to_title_easy <- function(pmid, date_type = c("PubDate", "pubmed", "entrez", "medline"), verbose = FALSE) {
#  pmid_to_title_easy <- function(pmid, date_type) {
# print(paste(pmid, date_type))
  # test: pmid <- "16571413"
  # test: pmid <- "21357945"
    # test: pmid <- "23844129"
    # test: pmid <- "22617845"
    # test: pmid <- "29868000"
    # test: pmid <- "28137280"
  pmids_thing <-  get_pubmed_ids(pmid)

  article <- fetch_pubmed_data(pmids_thing, format = "xml")

  #article <- fetch_pubmed_data(pmids_thing, format = "medline")
  # Get real publication year and month the hard way
  x <- read_xml(article)  # returns "xml_document" "xml_node"
  x <- iconv(x)  # returns character
  # At each newline, a new list item is added.  Example of a date entry (hours and minutes removed):
  # <PubMedPubDate PubStatus=\"pubmed\">\n <Year>2013</Year>\n <Month>5</Month>\n <Day>10</Day>\n </PubMedPubDate>\n
  x <- strsplit(x, split = "\n", fixed = TRUE)  # returns list

  # x is a list of length 1 of N character vectors
  # We can get rid of the top level list
  x <- x[[1]]
  x <- lapply(x, trimws)

  # returns matrix[n,1], change to vector, though not necessary
  s <- as.vector(sapply(x, function(y) {grepl("<AbstractText", y)}))
  w <- which(s)
  abstract <- sapply(x[w], function(y) {
    custom_grep(y, tag = "AbstractText", format = "character")
  })

  tokens <- sapply(x[w], function(y) {
    strsplit(y, split = " ")
  })

  li <- sapply(tokens, function(x) {
    p <- grep("Label=\"", x)
    if(length(p) == 0) {
      p <- NA
    }
    return(p)
  })

  labels <- mapply(function(x, n) {
    if(is.na(n)) {
      y <- "" # no label this section
    } else {
      y <- x[n]
      y <- sub("Label=\"", "", y)
      y <- sub("\\\".*", ":", y)
      if(substr(y, nchar(y), nchar(y)) != ":") {
        y <- paste0(y, ":")
      }
    }
    return(y)
  }, tokens, li)

  abstract <- mapply(function(x, n) {
    paste(x, n)
  }, labels, abstract)
  # this will also insert a blank space before a section with no label.
  abstract <- paste(abstract, collapse = " ")
  abstract <- trimws(abstract)

  if(!validUTF8(abstract)) {
    print(paste("PMID", pmid, "has non-UTF-8 character(s)"))
  }

    # returns matrix[n,1], change to vector, though not necessary
  # date_type must be enclosed in escaped quotes
  if (date_type == "PubDate") {
    s <- as.vector(sapply(x, function(x) {grepl("<PubDate>", x)}))
  } else {
    s <- as.vector(sapply(x, function(x) {grepl(paste0("PubStatus=\"", date_type, "\">"), x)}))
  }

  w <- which(s)
  custom_grep(x[w+1], tag = "Year", format = "character")
  pubYear <- custom_grep(x[w+1], tag = "Year", format = "character")
  if(is.null(pubYear)) {
    pubYear <- getPMIDYear(pmid)
    if(pubYear == "") {
      pubYear <- "0"
      if (verbose) {
        print(paste("No year found for PMID:", pmid))
      }
    }
  }

  pubMonth <- custom_grep(x[w+2], tag = "Month", format = "character")
  if(is.null(pubMonth)) {
    pubMonth <- getPMIDMonth(pmid)
    if(pubMonth == "") {
      pubMonth <- "0"
      if (verbose) {
        print(paste("No month found for PMID:", pmid))
      }
    }
  }
  # For PubDate, the month can be MMM, NN, or just missing...
  monthInt <- match(pubMonth, month.abb)
  if(!is.na(monthInt)) {
    pubMonth <- monthInt
  }
  if(pubMonth %in% 0:9) {
    pubMonth <- paste0("0", pubMonth)
  }

  pubDay <- custom_grep(x[w+3], tag = "Day", format = "character")
  if(is.null(pubDay)) {
    pubDay <- "00"
    if (verbose) {
      print(paste("No day found for PMID:", pmid))
    }
  } else if(pubDay %in% 0:9) {
    pubDay <- paste0("0", pubDay)
  }

  epm <- table_articles_byAuth(pubmed_data = article)

  # another way
  # epm <- article_to_df(article)  # does not return actual publication date

  if (length(epm$lastname) > 1)  {
    aname <- paste(epm$lastname[1], "et al.")  # each author's last name is in successive entries, just take first.
  } else {
    aname <- epm$lastname
  }

  auth_ref <- paste0(aname,  " (", pubYear, ")")  # replace epm$year[1] with pubYear
  pubmed_ref <- paste0("PMID: <a href = 'https://www.ncbi.nlm.nih.gov/pubmed/?term=", pmid, "' target='_blank'>", pmid, "</a>")

  title <- paste(auth_ref, epm$title[1], epm$jabbrv[1], pubmed_ref)

  # date <- paste(epm$year[1], epm$month[1], epm$day[1], sep = ".")
  if (date_type == "PubDate") {
    date <- paste(pubYear, pubMonth, sep = ".")
  } else {
    date <- paste(pubYear, pubMonth, pubDay, sep = ".")
  }
  return(data.frame(pmid, author = epm$lastname[1], title, date, abstract, stringsAsFactors = FALSE))
}

