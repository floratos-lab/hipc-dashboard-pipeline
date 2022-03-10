# It would be nice to use a direct call like this, then parse the xml return.
# url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=17651872"
# text <- read_html(url)
# But instead, try using easyPubMed.

# NOTE - There are a number of hard-coded special cases in the lookup tables below.

library(xml2)
library(easyPubMed)

### DATA ###
# Special cases for consortia as first author:
# Actual author names included for reference
# Look at the actual journal article to see what the first author was.
# If the consortium was the author, set "use" to TRUE.
use_consortium <- rbind(c("28854372", "Rechtien", FALSE),
                        c("23759749", "Richert", FALSE),
                        c("30244872", "Fischer", FALSE),
                        c("28404856", "Huttner", FALSE),
                        c("28842433", "HIPC-CHI Signatures Project Team", TRUE),
                        c("24725414", "Tsang", FALSE),
                        c("32826343", "Maucourant", FALSE),
                        c("32717743", "Lucas", FALSE),
                        c("33846275", "Saris", FALSE),
                        c("32856707", "Young", FALSE))

use_consortium <- as.data.frame(use_consortium, stringsAsFactors = FALSE)
colnames(use_consortium) <- c("pmid", "author", "use")
rownames(use_consortium) <- use_consortium$pmid
use_consortium$use <- as.logical(use_consortium$use)


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
# Dropped support for PubDate and medline.
# Will just use "pubmed" now.
pmid_to_title_easy <- function(pmid, print_pub_year) {
#  pmid_to_title_easy <- function(pmid, date_type) {
# print(paste(pmid, date_type))
  # test: pmid <- "16571413"
  # test: pmid <- "21357945"
    # test: pmid <- "23844129"
    # test: pmid <- "22617845"
    # test: pmid <- "29868000"
    # test: pmid <- "28137280"
    # test: pmid <- "19155521" # has name Garcia with special characters
    # test: pmid <- "33033248"
    # test: pmid <- "32826343"
  pmids_thing <-  get_pubmed_ids(pmid)

  article <- fetch_pubmed_data(pmids_thing, format = "xml")

  # Get real publication year and month the hard way
  x <- read_xml(article)  # returns "xml_document" "xml_node"
  # FIXME - "article" has the correct special UTF-8 characters for e.g. accented names.
  #         They are lost here.
  x <- iconv(x)  # returns character
  # At each newline, a new list item is added.  Example of a date entry (hours and minutes removed):
  # <PubMedPubDate PubStatus=\"pubmed\">\n <Year>2013</Year>\n <Month>5</Month>\n <Day>10</Day>\n </PubMedPubDate>\n
  x <- strsplit(x, split = "\n", fixed = TRUE)  # returns list

  # x is a list of length 1 of N character vectors
  # We can get rid of the top level list
  x <- x[[1]]
  x <- lapply(x, trimws)


    # returns matrix[n,1], change to vector, though not necessary
  # date_type must be enclosed in escaped quotes

  s <- as.vector(sapply(x, function(x) {grepl(paste0("PubStatus=\"", "pubmed", "\">"), x)}))

  w <- which(s)
  custom_grep(x[w+1], tag = "Year", format = "character")
  pubYear <- custom_grep(x[w+1], tag = "Year", format = "character")
  pubMonth <- custom_grep(x[w+2], tag = "Month", format = "character")

  if(pubMonth %in% 0:9) {
    pubMonth <- paste0("0", pubMonth)
  }

  pubDay <- custom_grep(x[w+3], tag = "Day", format = "character")
  if(pubDay %in% 0:9) {
    pubDay <- paste0("0", pubDay)
  }

  # Note - table_articles_byAuth() failed for article with PMID 28842433,
  # which has a collectiveName value. No problem with other such articles.
  s <- as.vector(sapply(x, function(y) {grepl("<CollectiveName>", y)}))
  w <- which(s)  # note that only first value in list x[w] used.
  if(length(w) > 0) {
    collective_name <- custom_grep(x[w], tag = "CollectiveName", format = "character")
    if(!(pmid %in% use_consortium$pmid)) {
      print(paste("New collective name found for PMID", pmid, collective_name))
    }
  } else {
    collective_name <- NULL
  }
  # if wanted all collective names
  #sapply(x[w], custom_grep, tag = "CollectiveName", format = "character")

  # Text retrieved using table_articles_byAuth() handles special characters properly.
  # But if it fails to return a value, then use direct parse method.
  epm <- table_articles_byAuth(pubmed_data = article)

  if(length(epm$lastname[1]) == 0) {
    print(paste("table_articles_byAuth() did not find author name; instead using direct parse"))

    s <- as.vector(sapply(x, function(y) {grepl("<LastName>", y)}))
    w <- which(s)[1]  # note that only first value in list x[w] used anyway.
    author_lname <- custom_grep(x[w], tag = "LastName", format = "character")
    if (length(s) > 1) {
      author_lname <- paste(author_lname, "et al.")  # each author's last name is in successive entries, just take first.
    }
  } else {
    author_lname <-  epm$lastname[1]
    if (length(epm$lastname) > 1) {
      author_lname <- paste(author_lname, "et al.")  # each author's last name is in successive entries, just take first.
    }
  }

  s <- as.vector(sapply(x, function(y) {grepl("<ISOAbbreviation>", y)}))
  w <- which(s)  # note that only first value in list x[w] used.
  journal_abb <- custom_grep(x[w], tag = "ISOAbbreviation", format = "character")

  if(length(epm$title) == 0) {
    print(paste("table_articles_byAuth() did not find title, using direct parse"))

    s <- as.vector(sapply(x, function(y) {grepl("<ArticleTitle>", y)}))
    w <- which(s)  # note that only first value in list x[w] used.
    article_title <- custom_grep(x[w], tag = "ArticleTitle", format = "character")
  } else {
    article_title <- epm$title[1]
  }

  aname <- author_lname # default
  if(pmid %in% use_consortium$pmid) {
    if(use_consortium[as.character(pmid), "use"] == TRUE) {
      if(!is.null(collective_name)  && length(collective_name) > 0) {
        aname <- collective_name
      }
    }
  }
  # Use externally supplied print_pub_year instead of pubmed pubYear value obtained from NCBI.
  auth_ref <- paste0(aname,  " (", print_pub_year, ")")  # replace epm$year[1] with pubYear
  pubmed_ref <- paste0("PMID: <a href = 'https://www.ncbi.nlm.nih.gov/pubmed/?term=", pmid, "' target='_blank'>", pmid, "</a>")
  dashboard_title <- paste(auth_ref, article_title, paste0(journal_abb, "."), pubmed_ref)

  date <- paste(pubYear, pubMonth, pubDay, sep = ".")
  abstract <- ""
  return(data.frame(pmid, author = aname, article_title, dashboard_title, date, abstract, stringsAsFactors = FALSE))
}

get_titles_and_dates <- function(df2, renew_pmids, pmid_file, log_files, pmid_logfile_basename) {
  pmids_uniq <- unique(df2$publication_reference_id)
  
  summary_df <- add_to_summary(summary_df, "Unique PMIDs", length(pmids_uniq))
  # Get user-entered year for each pmid
  # FIXME - need to catch if user did not enter year.
  print_pub_year <- df2$publication_date[match(pmids_uniq, df2$publication_reference_id)]
  # FIXME - need to check date format valid
  print_pub_year <- stri_sub(print_pub_year, 1, 4)
  
  if (renew_pmids || !file.exists(pmid_file)) {
    td <- mapply(pmid_to_title_easy, pmids_uniq, print_pub_year, SIMPLIFY = FALSE)
    titles_and_dates_df <- as.data.frame(rbindlist(td))
    save(titles_and_dates_df, file = pmid_file)
    write.table(titles_and_dates_df,
                file = paste(log_files,
                             paste(pmid_logfile_basename, "titles_and_dates_df.tsv", sep = "-"),
                             sep = "/"), 
                row.names = FALSE, col.names = TRUE, sep = "\t")
  } else {
    load(file = pmid_file)
  }
  return(titles_and_dates_df)
}

