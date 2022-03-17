#' Converts Idenifiers to pathway level for species translation. 
#'
#' Can take. Entrez, Ensembl gene, Ensembl transcript, Entrez, RefSeq and Gene symbol. Is using Entrez as key to the pathways
#' @import magrittr
#' @importFrom GSEABase getGmt
#' @importFrom singscore  multiScore
#' importFrom singscore, rankGenes)
#' importFrom dplyr select group_by_at across everything summarise filter one_of
#' @importFrom tibble column_to_rownames
#' @param counts An expression matrix of class matrix
#' @param species Options are "Human" (default), "Mouse", and "RAT".
#' @return An expression matrix pathway scores as row names.
#' @noRd
#' @examples
#' convertGeneName(exprs_dmap[1:20,1:10],to="hgnc_symbol")

convertToPathway <- function(counts, species){
        # remove version number
        cur_genes <- gsub("\\.[0-9]+$","",rownames(counts))

        # get gene mapper data frame
        gene_mapper <- mapGene(cur_genes,to="entrez")

        # case1: mapGene fails to infer gene id
        if (is.null(gene_mapper)){
                stop("Could not infer gene identifiers from row names of expression matrix.
                     Please set annotation = NA to stop gene id conversion and make sure manually that new samples have same gene ids as scaffold dataset")
        }

        from <- colnames(gene_mapper)[1]

        # case 2: from==to. Only need to remove version number before return
        if (from==to) {
                rownames(counts) <- cur_genes
                return(as.matrix(counts))
        }

        # case 3: conversion between gene id and gene id or transcript id. Add all matches
        # add all matches
        counts <- merge(gene_mapper,counts,by.x=from,by.y=0,all=F)
        counts <- counts[,-1]
        counts <- counts%>% dplyr::group_by_at(to) %>%
                dplyr::summarise(dplyr::across(dplyr::everything(),sum)) %>%
                tibble::column_to_rownames(var=to)
        

        gene_sets <- getGmt("data/c5.go.bp.v7.3.entrez.gmt")
        scores <- multiScore(rankData = rankGenes(expr), upSetColc = gene_sets, knownDirection = TRUE)
        scores$Scores

        return (as.matrix(scores$Scores))

}
