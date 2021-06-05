#' Create ExpressionSet object
#'
#' This function creates ExpressionSet from an expression matrix and a phenotype data frame.
#' This function also performs data preprocessing: samples lacking either count data or phenotype annotation are removed. NAs are replaced by 0.
#'
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame
#' @param counts An expression matrix of class matrix, or data frame that can be converted to matrix.
#' @param pheno A dataframe corresponding to the expression matrix.
#' @param colname A column name of \code{pheno}. Differential analysis will be performed using this column of phenotype as independent variables.
#' @return An ExpressionSet object
#' @noRd
#' @examples
#' create_eset(exprs_dmap,pData_dmap)

createEset <- function(counts,pheno,colname,to="ensembl_gene_id"){

        # check class of counts, convert to matrix
        if (!is.matrix(counts)){
                tryCatch(counts<- data.matrix(counts),
                         error=function(e) {
                                 stop("Expression matrix cannot be converted to matrix:",e)
                         })
        }

        # check class of pheno, convert to dataframe
        if (!is.data.frame(pheno)){
                tryCatch(pheno <- as.data.frame(pheno),
                         error=function(e){
                                 stop("Phenotype data cannot be converted to data frame:", e)
                         })
        }

        # convert gene names
        counts <- convertGeneName(counts,to=to)

        #select the specified column of phenotype table as final phenotype table
        pheno <- pheno[,colname,drop=F]

        # remove rows with NA in pheno and throw message
        complete_idx <- which(complete.cases(pheno))
        na_rownum <- dim(pheno)[1]-length(complete_idx)
        if (na_rownum){
                pheno <- pheno[complete_idx,,drop=F]
                message(na_rownum, " row(s) in phenotype table have no data in the required column, thus removed.")
        }

        # subset counts to only contain samples with annotation
        idx <- which(!colnames(counts) %in% rownames(pheno))
        if (length(idx)>0){
                counts <- counts[,-idx,drop=F]
                message(length(idx)," sample(s) have no phenotype annotation, thus removed.")

        }
        idx <- which(! rownames(pheno) %in% colnames(counts))
        if (length(idx)>0){
                pheno <- pheno[,-idx,drop=F]
                message(length(idx)," sample(s) have no expression count data, thus removed.")

        }

        # match row/col names
        idx <- match(colnames(counts),rownames(pheno))
        pheno <- pheno[idx,,drop=F]

        # fill NA in count matrix with 0.
        sum_na <-sum(is.na(counts))
        if (sum_na>0){
                counts[is.na(counts)] <- 0
                message("Count matrix has ",sum_na, " missing values. Replace NA by 0.")
        }

        # return ExpressionSet
        Biobase::ExpressionSet(counts, phenoData=Biobase::AnnotatedDataFrame(pheno))


}
