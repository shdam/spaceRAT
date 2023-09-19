#' Check if input is a prebuilt scaffold
#'
#' @noRd
checkPrebuilt <- function(object, classes){
    stopifnot("Please only check a character string" = is(object, "character"))
    if( # Return prebuilt scaffold
        object %in% listScaffolds() && is(classes, "NULL")
        ){
        space <- loadData(object)
        return(space)
    } else if( # Subsetting a Scaffold
        object == "DMAP_scaffold" && !is(classes, "NULL")
        ){
        object <- loadData("DMAP_exprs")
        pheno <- loadData("DMAP_pData")
        colname <- "cell_types"
    } else{ # Unknown scaffold
        stop(
            "Incorrectly named prebuilt scaffold. The available are: ",
            listScaffolds()
        )
    }
    # Return subsetted scaffold
    return(list("object" = object, "pheno" = pheno, "colname" = colname))
}
