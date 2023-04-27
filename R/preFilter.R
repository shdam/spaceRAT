#' Prefilter counts data
#'
#' @inheritParams buildScaffold
#' @noRd
preFilter <- function(counts, data = "logged", threshold = 10){
    # data preprocessing and create eset
    if (data == "logged"){
        # remove genes with total count<10
        idx <- which(rowSums(exp(counts))<threshold)
        if (length(idx)==dim(counts)[1])
            stop("Low quality data!All genes have total counts less than 10.")
        if (length(idx)>0) counts <- counts[-idx,]

    } else if(data == "raw"){
        # ensure no negative value
        if (any(counts<0))
            stop("Negative values are not allowed in raw count matrix!")
        # remove genes with total count<10
        idx <- which(rowSums(counts)<threshold)
        if (length(idx)==dim(counts)[1])
            stop("Low quality data! All genes have total counts less than 10.")
        if (length(idx)>0) counts <- counts[-idx,]
        counts <- log(counts+1)
    }
    return(counts)
}
