#' Convert between Ensembl gene, Ensembl transcript, Entrez, RefSeq and Gene symbol for character vector of gene identifiers
#'
#' This function takes in a character vector of gene names, determine which gene identifier it is, then converts the gene names to the specified type, finally returns data frame as gene mapper.
#' @import magrittr
#' @importFrom dplyr select group_by_at across everything summarise filter
#' @param vec Character vector of gene identifier.
#' @param to A character specifying the gene identifier to convert to. Options are "ensembl_gene", "ensembl_transcript", "entrez", "hgnc_symbol" and "refseq_mrna".
#' @return A data frame with original gene identifier as first column, converted gene identifier as second column.
#' @noRd

mapGene <- function(vec,to){

        vec <- as.character(vec)
        from <- NULL

        # determine which gene id is used
        if(all(startsWith(vec,"ENSG"))){
                from <- "ensembl_gene"
        }else if(all(startsWith(vec,"ENST"))){
                from <- "ensembl_transcript"
        }else if(all(startsWith(vec,"NM_"))){
                from <- "refseq_mrna"
        }else if(all(grepl("^[0-9]+$",vec))){
                from <- "entrez"
        }else{
                idx <- which(gene_id_converter_hs$hgnc_symbol %in% vec)
                if (length(idx)>0) from <- "hgnc_symbol"
        }

        # fail to determine gene id
        if (is.null(from)) {
                return(NULL)

        }else if (from==to) {
                df <- data.frame(from=vec,to=vec)
                colnames(df) <- c(from,to)
                return(df)

        }else if(from!=to){
                message("Convert gene indentifiers of count matrix from ",from," to ",to,".")
        }

        # subset gene_id_converter_hs
        # if from=="hgnc_symbol", then the idx has already be calculated. Avoid repetitive calculation.
        if(from!= "hgnc_symbol") {
                idx <- which(gene_id_converter_hs[[from]] %in% vec)
        }

        df <- gene_id_converter_hs[idx,c(from,to),drop=F]
        df <- df[!duplicated(df),]
        df <- df[complete.cases(df),]
        df <- df %>% dplyr::group_by_at(from) %>%
                dplyr::filter(row_number()==1)%>%
                as.data.frame()

        return(df)
}
