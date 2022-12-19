#' Check matrix has correct format
#'
#'
#' @param mat A matrix to check
#' @noRd
checkMatrix <- function(mat){

  # Identify non-numerics in counts data
  numerics <- vapply(mat[seq_len(5),], is, "numeric", FUN.VALUE = TRUE)
  stopifnot("Too many non-numerics in counts data! Please separate from counts." = sum(!numerics) < 2)
  empty_rownames <- all(rownames(mat)[seq_len(6)] == paste(seq_len(6)))

  # No rowname information
  stopifnot("The expression matrix does not contain rownames" = (sum(!numerics) == 1 || !empty_rownames))

  # Check counts rownames
  if(sum(!numerics) == 1){

    if(empty_rownames){ # No rownames
      if(is(mat, "data.frame")){
        mat <- tibble::column_to_rownames(mat, colnames(mat)[which(!numerics)])
      } else if(is(mat, "matrix")){
        rownames(mat) <- mat[which(!numerics)]
        mat <- mat[which(numerics)]
      }
    } else if(!empty_rownames){ # Rownames and additional column
      warning("column: ", colnames(mat)[which(!numerics)], " is removed from expression matrix. If that column was intended to be the matrix' rownames, please manually set the rownames.")
      mat <- mat[which(numerics)]
    }
  } # Else the matrix is fine

  return(as.matrix(mat))

}
