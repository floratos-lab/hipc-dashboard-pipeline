# It would be nice to use a direct call like this, then parse the xml return.
# url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=17651872"
# text <- read_html(url)
# But instead, try using easyPubMed.

# NOTE - There are a number of hard-coded special cases in the lookup tables below.

library(xml2)
library(easyPubMed)
library(stringi)

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
  pmids_thing <-  get_pubmed_ids(pmid)

  article <- fetch_pubmed_data(pmids_thing, format = "xml")

  # Get real publication year and month the hard way
  x <- read_xml(article)  # returns "xml_document" "xml_node"
  pub_dates <- xml_find_all(x, ".//PubMedPubDate")
  pubmed_date <- pub_dates[xml_attr(pub_dates, "PubStatus")=="pubmed"]
  pubYear <- xml_text(xml_child(pubmed_date, "Year"))
  pubMonth <- xml_text(xml_child(pubmed_date, "Month"))
  pubDay <- xml_text(xml_child(pubmed_date, "Day"))

  if(pubMonth %in% 0:9) {
    pubMonth <- paste0("0", pubMonth)
  }

  if(pubDay %in% 0:9) {
    pubDay <- paste0("0", pubDay)
  }

  collective_name <- xml_text(xml_find_first(x, ".//CollectiveName"))
  lastnames <- xml_find_all(x, ".//LastName")
  author_lname <- paste0(xml_text(lastnames[[1]]), ifelse(length(lastnames)>1, " et al.", ""))
  journal_abb <- xml_text(xml_find_first(x, ".//ISOAbbreviation"))
  article_title <- xml_text(xml_find_first(x, ".//ArticleTitle"))

  aname <- author_lname # default
  if(pmid %in% use_consortium$pmid) {
    if(use_consortium[as.character(pmid), "use"] == TRUE) {
      if(!is.null(collective_name)  && length(collective_name) > 0) {
        aname <- collective_name
      }
    }
  }
  # Use externally supplied print_pub_year instead of pubmed pubYear value obtained from NCBI.
  auth_ref <- paste0(aname,  " (", print_pub_year, ")")
  pubmed_ref <- paste0("PMID: <a href = 'https://www.ncbi.nlm.nih.gov/pubmed/?term=", pmid, "' target='_blank'>", pmid, "</a>")
  dashboard_title <- paste(auth_ref, article_title, paste0(journal_abb, "."), pubmed_ref)

  date <- paste(pubYear, pubMonth, pubDay, sep = ".")
  return(data.frame(pmid, dashboard_title, date, stringsAsFactors = FALSE))
}
