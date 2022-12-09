#' Check matrix has correct format
#'
#'
#'
checkMatrix <- function(mat){

  # Identify non-numerics in counts data
  numerics <- sapply(mat[1:5,], is, "numeric")
  stopifnot("Too many non-numerics in counts data! Please separate from counts." = sum(!numerics) < 2)
  empty_rownames <- all(rownames(mat)[1:6] == paste(1:6))

  # No rowname information
  stopifnot("The expression matrix does not contain rownames" = (sum(!numerics) == 1 || !empty_rownames))

  # Check counts rownames
  if(sum(!numerics) == 1){

    if(empty_rownames){ # No rownames
      if(is(mat, "data.frame")){
        mat <- tibble::column_to_rownames(mat, colnames(mat)[which(!numerics)])
      } else if(is(mat, "matrix")){
        rownames(mat) <- counts[which(!numerics)]
        mat <- mat[which(numerics)]
      }
    } else if(!empty_rownames){ # Rownames and additional column
      warning(paste("column:", colnames(mat)[which(!numerics)], "is removed from expression matrix. If that column was intended to be the matrix' rownames, please manually set the rownames."))
      mat <- mat[which(numerics)]
    }
  } # Else the matrix is fine

  return(as.matrix(mat))

}
