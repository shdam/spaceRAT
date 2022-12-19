#' Convert between Ensembl gene, Ensembl transcript, Entrez, RefSeq,
#' and Gene symbol for count ?matrix
#'
#' This function takes in a count matrix, and converts its row names to
#' Ensembl (by default), or any specified gene identifier.
#' @param counts An expression matrix of class matrix
#' @param to A character specifying the gene id to convert to.
#' Options are "ensembl_gene" (default), "ensembl_transcript", "entrez",
#' "hgnc_symbol", and "refseq_mrna".
#' @return A count matrix with transformed gene names as row names.
#' @noRd
#' @examples
#' utils::data("exprs_dmap")
#' convertGeneName(exprs_dmap[seq_len(20),seq_len(10)],to="hgnc_symbol")

convertGeneName <- function(counts, to = "ensembl_gene"){

# remove version number
cur_genes <- gsub("\\.[0-9]+$","", rownames(counts))

    # get gene mapper data frame
    gene_mapper <- mapGene(cur_genes, to = to)

    # case1: mapGene fails to infer gene id
    if (is.null(gene_mapper)){
        stop("Could not infer gene identifiers from row names
        of expression matrix.
        Please set annotation = NULL to stop gene id conversion and make sure
             manually that new samples have same gene ids as scaffold dataset")
        }

    from <- colnames(gene_mapper)[1]

    # case 2: from==to. Only need to remove version number before return
    if (from == to) {
        rownames(counts) <- cur_genes
        return(counts)
    }

    # case 3: conversion between gene id and gene id or transcript id.
    # Add all matches
    matches <- match(cur_genes, gene_mapper[[from]], nomatch = 0)
    cur_genes[which(matches != 0)] <- gene_mapper[[to]][matches]

    if(any(matches == 0)){
        counts <- counts[which(matches != 0), ]
        cur_genes <- cur_genes[which(matches != 0)]
    }

    if(length(unique(cur_genes)) != length(cur_genes)){
        # Sum genes with the same new id
        extra_warning <- "\nOBS: Duplicate gene names were found after
        conversion. Expression values will be summed."

        # Sum duplicates
        counts <- stats::aggregate(counts, by = list(cur_genes),
                                   FUN = sum) %>%
            tibble::column_to_rownames(var = "Group.1")

    } else{
        extra_warning <- ""
        rownames(counts) <- cur_genes
    }

    if(any(matches == 0)) warning(
        "Be aware only ",
        round(sum(matches > 0)/length(matches)*100, 2),
        "% of gene names were converted from ", from, " to ", to,
        ". The non-converted genes will be removed.
        Use 'annotation = NULL' to prevent conversion.",
        extra_warning
    )

    # counts <- merge(gene_mapper, counts, by.x = from, by.y = 0,
    #                 all = FALSE)
    # counts <- counts[, -1]
    # counts <- counts %>% dplyr::group_by(.data[[to]]) %>%
    #         dplyr::summarise(dplyr::across(dplyr::everything(),sum)) %>%
    #         tibble::column_to_rownames(var = to)
    return (counts)
}
