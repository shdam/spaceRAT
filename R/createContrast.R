#' Create contrast matrix
#'
#' This function takes in a ExpressionSet object, count the number of each cell type, then create a contrast matrix to contrast each cell type (with at least 3 cells) against all other cell types.
#' @importFrom Biobase pData
#' @param eset an ExpressionSet
#' @export
#' @return A contrast matrix
#' @examples
#' eset_dmap <- createEset(exprs_dmap,pData_dmap,"cell_types")
#' createContrast(eset_dmap)
createContrast <- function(eset){
        count_table <- data.frame(table(Biobase::pData(eset)[,1]))
        cell_types <- count_table[count_table$Freq>2,]$Var1 # create contrast only for cell types with at least 3 cells
        n <- length(cell_types)
        mat <- matrix(-1/(n-1),n,n)
        diag(mat) <- 1
        rownames(mat) <- cell_types
        colnames(mat) <- paste(cell_types,"others",sep="_")
        mat
}
