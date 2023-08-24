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
#' @importFrom stats aggregate
#' @noRd
convertGeneName <- function(object, assay = "counts", to = "ensembl_gene") {
    # Ensure the provided object is a SummarizedExperiment
    if (!is(object, "SummarizedExperiment")) {
        stop("Please provide a SummarizedExperiment object.")
    }

    # remove version number
    cur_genes <- gsub("\\.[0-9]+$","", rownames(object))

    # get gene mapper data frame
    gene_mapper <- mapGene(cur_genes, to = to)

    # case1: mapGene fails to infer gene id
    if (is.null(gene_mapper)){
        stop("Could not infer gene identifiers from rownames ",
             "of expression matrix. Please set annotation = NULL to stop ",
             "gene id conversion and make sure manually that new samples ",
             "have same gene ids as scaffold dataset")
    }

    from <- colnames(gene_mapper)[1]
    # case 2: from==to. Only need to remove version number before return
    if (from == to) {
        rownames(object) <- cur_genes
        return(object)
    }

    # case 3: conversion between gene id and gene id or transcript id.
    matches <- match(cur_genes, gene_mapper[[from]], nomatch = 0)
    cur_genes[which(matches != 0)] <- gene_mapper[[to]][matches]


    if(any(matches == 0)){
        object <- object[which(matches != 0), ]
        cur_genes <- cur_genes[which(matches != 0)]
    }

    if(length(unique(cur_genes)) != length(cur_genes)){
        # Access the counts from the SummarizedExperiment
        counts <- assay(object, assay = assay)

        # Sum genes with the same new id
        extra_warning <- paste(
            "\nOBS: Duplicate gene names were found after ",
            "conversion. Expression values will be summed.")
        counts <- stats::aggregate(counts, by = list(cur_genes), FUN = sum)
        rownames(counts) <- counts$Group.1
        counts$Group.1 <- NULL

        # Set the updated counts back to the SummarizedExperiment object
        metadata <- metadata(object)
        coldata <- colData(object)
        object <- SummarizedExperiment(
            assays = list("counts" = counts),
            colData = coldata,
            metadata = metadata
        )
        assay(object, assay) <- counts
    } else {
        extra_warning <- ""
        rownames(object) <- cur_genes
    }


    if(any(matches == 0)) warning(
        "Be aware only ",
        round(sum(matches > 0)/length(matches)*100, 2),
        "% of gene names were converted from ", from, " to ", to,
        ". The non-converted genes will be removed.
        Use 'annotation = NULL' to prevent conversion.",
        extra_warning)

    return(object)
}
