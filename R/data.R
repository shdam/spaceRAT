#' @title Expression matrix for scaffold
#' @source Kindly offered by Fred
#' @description  exprs_map records expression matrix of normal hematopoietic cells. The data is measured by microarray.
#' @docType data
#' @usage exprs_dmap
#' @format a matrix
#' @keywords datasets
"exprs_dmap"

#' @title Phenotype data for scaffold
#' @source Kindly offered by Fred
#' @description  pData_map records phenotype information corresponding to exprs_dmap.
#' @docType data
#' @usage pData_dmap
#' @format a data frame
#' @keywords datasets
"pData_dmap"


#' @title Expression matrix as new samples
#' @source Kindly offered by Fred
#' @description exprs_ilaria contains expression matrix of Acute Myeloid Leukemia of erythroid subtype.
#' @docType data
#' @usage exprs_ilaria
#' @format a data frame
#' @keywords datasets
"exprs_ilaria"


#' @title Phenotype data for new samples
#' @source Kindly offered by Fred
#' @description pData_ilaria contains phenotype data corresponding to exprs_ilaria
#' @docType data
#' @usage pData_ilaria
#' @format a data frame
#' @keywords datasets
"pData_ilaria"

#' @title  Gene name mapper for human
#' @source downloaded from biomaRt
#' @description This table contains ensembl_gene_id of all human genes, and corresponding ensembl_transcript_id, entrezgene_id ,hgnc_symbol, and refseq_mrna
#' @docType data
#' @usage gene_id_converter_hs
#' @format a data frame
#' @keywords datasets
"gene_id_converter_hs"

#' @title  Gene name mapper for Mus musculus (mouse)
#' @source downloaded from biomaRt
#' @description This table contains ensembl_gene_id of all mouse genes, and corresponding ensembl_transcript_id, entrezgene_id, and mgi_symbol
#' @docType data
#' @usage gene_id_converter_mm
#' @format a data frame
#' @keywords datasets
"gene_name_mapper_mm"

#' @title  DMAP scaffold
#' @source Created from exprs_dmap
#' @description This scaffold was build from the \code{\link{exprs_dmap}} object.
#' @docType data
#' @usage DMAP_scaffold
#' @format a scaffoldSpace
#' @keywords datasets
"DMAP_scaffold"
