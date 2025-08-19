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
convertGeneName <- function(mat, to = "ensembl_gene") {

    if (grepl("^NM", rownames(mat)[1])) {
        cur_genes <- rownames(mat)
    } else {
        cur_genes <-  gsub("_.*$", "", rownames(mat))
    }

    # if (!grepl("^NM", rownames(mat)[1])) cur_genes <- gsub("_.*$", "", rownames(mat)) else cur_genes <- rownames(mat)
    cur_genes <- gsub("\\.[0-9]+$","", cur_genes) # remove version number

    gene_mapper <- mapGene(cur_genes, to = to) # get gene mapper data frame

    if (is.null(gene_mapper)){ # case1: mapGene fails to infer gene id
        stop(
            "Could not infer gene identifiers from rownames ",
            "of expression matrix. Please set annotation = NULL to stop ",
            "gene id conversion and make sure manually that new samples ",
            "have same gene ids as scaffold dataset")
    }
    from <- colnames(gene_mapper)[1]
    # case 2: from==to. Only need to remove version number before return
    if (from == to) {
        rownames(mat) <- cur_genes
        return(mat)
    }
    # case 3: conversion between gene id and gene id or transcript id.
    matches <- match(cur_genes, gene_mapper[[from]], nomatch = 0)
    cur_genes[which(matches != 0)] <- gene_mapper[[to]][matches]

    if (grepl("_.*$", rownames(mat)[1]) & !grepl("^NM", rownames(mat)[1])) cur_genes <- paste0(
        cur_genes, sub("^.*?(_.*$)", "\\1", rownames(mat)))

    if(any(matches == 0)){ # Remove nomatches
        mat <- mat[which(matches != 0), ]
        cur_genes <- cur_genes[which(matches != 0)]
    }
    if(length(unique(cur_genes)) != length(cur_genes)){
        # Sum genes with the same new id
        extra_warning <- paste(
            "\nOBS: Duplicate gene names were found after ",
            "conversion. Expression values will be summed.")
        mat <- stats::aggregate(mat, by = list(cur_genes), FUN = sum)
        rownames(mat) <- mat$Group.1
        mat$Group.1 <- NULL
    } else {
        extra_warning <- ""
        rownames(mat) <- cur_genes
    }
    if(any(matches == 0)) warning(
        "Be aware only ",
        round(sum(matches > 0)/length(matches)*100, 2),
        "% of gene names were converted from ", from, " to ", to,
        ". The non-converted genes will be removed.
        Use 'annotation = NULL' to prevent conversion.",
        extra_warning)
    return(mat)
}
