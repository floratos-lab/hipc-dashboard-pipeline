# this script is based on copying a portion of main_inf_gene.R
# for download PubMed info only.
# Thus all other irrelevant code is removed.

source("pmid_to_title_easy.R") # for pmid_to_title_easy

# Assume executing from the ./code directory
source_curations       <- "../data/source_curations"
reference_files        <- "../data/reference_files"

exclude_pmid <- "33361761"  # PMID(s) that are not suitable for processing

##### Set up file and template name components #####
first_data_row <- 9  # not counting the first row, which becomes a column header.

sheet_file     <- "hipc_vaccine - cell_type_frequency.tsv"

# other than these two lines, this script is the same as the one for vac_gene
base_filename  <- "vac_cell_type"
response_behavior_type_var <- "cell type frequency"

# only three columns are used: publication_reference_id, publication_date, response_behavior_type
relevant_columns <- c("publication_reference_id","publication_date", "response_behavior_type")

insub <- read.delim(file =  paste(source_curations, sheet_file, sep = "/"),
                    strip.white = TRUE,
                    stringsAsFactors = FALSE)
insub  <- insub[names(insub) %in% relevant_columns]

df2 <- insub[first_data_row:nrow(insub),]

# The "pmid:" tag is not currently required, but check if anything was entered
df2$publication_reference_id <- sub("pmid:", "", df2$publication_reference_id)
s <- sapply(as.integer(df2$publication_reference_id), is.integer)
if (!all(s)) {
  stop("non-numeric pmid")
}

# Infection templates contain more than one response_behavior_type, filter out all but current type
message("number of row before filtering by response type: ", nrow(df2))
df2 <- df2[df2$response_behavior_type == response_behavior_type_var, ]

# remove signatures from unusable PMIDs
df2 <- df2[!(df2$publication_reference_id %in% exclude_pmid), ]
message("number of rows used: ", nrow(df2))

# only keep the unique rows
unique_pmid <- unique(df2$publication_reference_id)
df2 <- df2[match(unique_pmid, df2$publication_reference_id), ]
message("number of unique pmids: ", nrow(df2))

print_pub_year <- substring(df2$publication_date, 1, 4)

td <- mapply(pmid_to_title_easy, df2$publication_reference_id, print_pub_year) # matrix of lists
attr_names <- rownames(td)
td <- matrix(unlist(td), nrow=nrow(td))
rownames(td) <- attr_names
titles_and_dates_df <- as.data.frame(t(td))

# there are six columns in this data frame: "pmid", "author", "article_title", "dashboard_title", "date", "abstract"
# only three are used: pmid, dashboard_title, date
titles_and_dates_df$author <- NULL
titles_and_dates_df$article_title <- NULL
titles_and_dates_df$abstract <- NULL

message("number of titles/dates: ", nrow(titles_and_dates_df))
save(titles_and_dates_df, file = paste0(reference_files, "/", base_filename, "-titles_and_dates_df.RData"))
