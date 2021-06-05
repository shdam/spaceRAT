#' Convert between Ensembl gene, Ensembl transcript, Entrez, RefSeq and Gene symbol for matrix
#'
#' This function takes in a count matrix, and convert its rownames to Ensembl (by default), or any specified gene name.
#' @import magrittr
#' @importFrom dplyr select group_by_at across everything summarise filter one_of
#' @importFrom tibble column_to_rownames
#' @param counts an expression matrix of class matrix
#' @param to a character specifying the gene id to convert to. Options are "ensembl_gene_id" (default), "ensembl_transcript_id", "entrezgene_id", "hgnc_symbol" and "refseq_mrna".
#' @return A count matrix with transformed gene names as row names
#' @noRd
#' @examples
#' convertGeneName(exprs_dmap,to="hgnc_symbol")

convertGeneName <- function(counts,to="ensembl_gene_id"){
        # remove version number
        cur_genes <- gsub("\\.[0-9]+$","",rownames(counts))
        gene_mapper <- mapGene(cur_genes,to=to)
        from <- colnames(gene_mapper)[1]
        if (from==to) {
                rownames(counts) <- cur_genes
                return(as.matrix(counts))
        }

        if (from %in% c("ensembl_gene_id","entrezgene_id","hgnc_symbol")){
                gene_mapper <- gene_mapper %>%
                      dplyr::group_by_at(to)%>%
                      dplyr::filter(row_number()==1) %>%
                      tibble::column_to_rownames(var=from)
                counts <- counts[rownames(gene_mapper),]
                rownames(counts) <- gene_mapper[,1]
                return(as.matrix(counts))
        }

        if (from %in% c("ensembl_transcript_id","refseq_mrna")){ # add all matches
                counts <- merge(gene_mapper,counts,by.x=from,by.y=0,all=F) %>%
                        dplyr::select(-dplyr::one_of(from))%>%
                        dplyr::group_by_at(to) %>%
                        dplyr::summarise(dplyr::across(dplyr::everything(),sum)) %>%
                        tibble::column_to_rownames(var=to)
                return (as.matrix(counts))
        }
}
