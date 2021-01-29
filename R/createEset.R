#' Create ExpressionSet object
#'
#' This function creates ExpressionSet from an expression matrix and a phenotype data frame
#'
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame
#' @param expr_mat an expression matrix of class matrix, or data frame that can be converted to matrix.
#' @param pheno a dataframe corresponding to the expression matrix. The row names of pheno must be identical to column names of expr_mat.
#' @export
#' @return An ExpressionSet object
#' @examples
#' create_eset(exprs_dmap,pData_dmap)

createEset <- function(expr_mat,pheno){

        # check class of exprs, convert to matrix
        if (!is.matrix(expr_mat)){
                tryCatch(expr_mat <- data.matrix(expr_mat),
                         error=function(e) {
                                 stop("Expression matrix cannot be converted to matrix:",e)
                         })
        }

        # check class of pData, convert to dataframe
        if (!is.data.frame(pheno)){
                tryCatch(pheno <- data.frame(pheno),
                         error=function(e){
                                 stop("Phenotype data cannot be converted to data frame:", e)
                         })
        }

        # check whether row/col names match
        if (any(! colnames(expr_mat) %in% rownames(pheno))){
                stop("Row/column names do not match! Some column name(s) of expression matrix doesn't exist in row names of phenotype data.")
        }
        if (any(! rownames(pheno) %in% colnames(expr_mat))) {
                stop("Row/column names do not match! Some row name(s) of phenotype data doesn't exist in column names of expression matrix.")

        }
        idx <- match(colnames(expr_mat),rownames(pheno))

        # match row/col names
        pheno <- pheno[idx,,drop=F]

        # return ExpressionSet
        Biobase::ExpressionSet(expr_mat, phenoData=Biobase::AnnotatedDataFrame(pheno))
}
