# this file is adopted from the original generate_HIPC_submissions.R

library(R.utils) # for gzip

save_standardized_curations <- function(df2, base_filename) {
    del_cols <- c(
        "submission_name", "submission_date",
        "template_name", "short_comment", "process_note"
    )
    df2tmp <- df2[!colnames(df2) %in% del_cols]
    df2tmp <- df2tmp[-1]

    filename <- paste0(
        "../data/standardized_curations/", base_filename,
        "-standardized_denormalized.tsv"
    )
    write.table(df2tmp,
        file = filename, sep = "\t",
        row.names = FALSE, col.names = TRUE
    )
    gzip(filename,
        destname = paste0(filename, ".gz"), overwrite = TRUE,
        remove = TRUE
    )
}

save_convenience_files <- function(
    df2, header_rows, base_filename,
    exposure_type, response_type) {
    if (exposure_type != "VACCINE" && exposure_type != "INFECTION") {
        stop("Incorrect exposure type encountered")
    }
    if (response_type != "GENE" && response_type != "CELLTYPE_FREQUENCY") {
        stop("Incorrect response type encountered")
    }

    if (response_type == "GENE") {
        response_behavior_type_var <- "gene expression"
    } else if (response_type == "CELLTYPE_FREQUENCY") {
        response_behavior_type_var <- "cell-type frequency"
    }

    convenience_files <- "../data/convenience_files/"

    uniq_sig_row_ids <- unique(df2$sig_row_id)
    resp_components_annotated <- vector("list", length(uniq_sig_row_ids))
    recreated_template <- vector("list", length(uniq_sig_row_ids))

    for (i in seq_along(uniq_sig_row_ids)) {
        df2tmp <- df2[df2$sig_row_id == uniq_sig_row_ids[i], ]
        # Recreate a full signature in one row
        base_row <- df2tmp[1, ] # get first row for this uniqID

        response_rowname <- paste(base_row$publication_reference_id,
            base_row$sig_subm_id, uniq_sig_row_ids[i],
            sep = "_"
        )
        response_description <- paste("PMID", base_row$publication_reference_id,
            response_behavior_type_var, base_row$sig_subm_id,
            sep = " "
        )

        # Use the full original set of response components
        # rather than just those for which a valid symbol was found.
        base_row$response_component_original <- paste(
            unique(df2tmp$response_component_original),
            collapse = "; "
        )

        base_row$exposure_material_id <- paste(
            unique(df2tmp$exposure_material_id),
            collapse = "; "
        )
        base_row$tissue_type_term_id <- paste(
            unique(df2tmp$tissue_type_term_id),
            collapse = "; "
        )

        if (response_type == "GENE") {
            base_row$response_component <- paste(
                unique(df2tmp$response_component),
                collapse = "; "
            )
            resp_components_annotated[[i]] <- c(
                response_rowname,
                response_description, unique(df2tmp$response_component)
            )
        } else if (response_type == "CELLTYPE_FREQUENCY") {
            full_sig <- unique(df2tmp$fully_qualified_response_component)
            # FIXME - only response_component is getting put back together?
            base_row$response_component <- paste(full_sig, collapse = "; ")
            base_row$response_component_id <- paste(
                unique(df2tmp$response_component_id),
                collapse = "; "
            )
            base_row$proterm_and_extra <- paste(
                unique(df2tmp$proterm_and_extra),
                collapse = "; "
            )
            base_row$fully_qualified_response_component <- paste(
                unique(df2tmp$fully_qualified_response_component),
                collapse = "; "
            )
            # The pro_ontology_id values are already separated by semicolons,
            # so change to commas
            # before potentially joining two lists of pro-terms.
            df2tmp$pro_ontology_id <- sapply(
                df2tmp$pro_ontology_id,
                function(x) {
                    gsub(";", ",", x)
                }
            )
            base_row$pro_ontology_id <- paste(
                unique(df2tmp$pro_ontology_id),
                collapse = "; "
            )

            resp_components_annotated[[i]] <- c(
                response_rowname, response_description, full_sig
            )
        }

        # Reconstitute target_pathogen and exposure_material_id
        if (exposure_type == "VACCINE") {
            base_row$target_pathogen_taxonid <- paste(
                unique(df2tmp$target_pathogen_taxonid),
                collapse = "; "
            )
        }

        recreated_template[[i]] <- base_row
    }

    names(resp_components_annotated) <- uniq_sig_row_ids

    # consolidate to a single data.frame
    recreated_template_df <- as.data.frame(rbindlist(recreated_template))
    if (any(colnames(header_rows) != colnames(recreated_template_df))) {
        stop("mismatch between header rows and recreated_template_df rows")
    }

    recreated_template_df <- rbind(header_rows, recreated_template_df)

    # First save a complete version for use in debugging/logging
    del_cols <- c("submission_name", "submission_date", "template_name")
    recreated_template_df <- recreated_template_df[
        !colnames(recreated_template_df) %in% del_cols
    ]

    # Set that first column name back to blank
    colnames(recreated_template_df)[1] <- ""

    del_cols <- c("sig_subm_id", "sig_row_id")

    recreated_template_df <- recreated_template_df[
        !colnames(recreated_template_df) %in% del_cols
    ]
    write.table(recreated_template_df,
        file = paste0(
            convenience_files,
            base_filename, "-standardized_curation_template.tsv"
        ),
        sep = "\t", row.names = FALSE
    )

    gmt_file <- paste0(
        convenience_files,
        base_filename, "-response_components.gmt.txt"
    )
    if (file.exists(gmt_file)) file.remove(gmt_file)
    lapply(
        resp_components_annotated,
        function(x) {
            write.table(paste(x, collapse = "\t"),
                file = gmt_file, row.names = FALSE, col.names = FALSE,
                quote = FALSE, append = TRUE
            )
        }
    )
    message("Finished creating convenience files")
}
