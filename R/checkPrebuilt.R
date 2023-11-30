#' Check if input is a prebuilt scaffold
#'
#' @noRd
checkPrebuilt <- function(object, classes){
    stopifnot("Please only check a character string" = is(object, "character"))
    if( # Return prebuilt scaffold
        object == "DMAP.v1" && !is(classes, "NULL")
        ){
        object <- loadData("DMAP_exprs")
        pheno <- loadData("DMAP_pData")
        colname <- "cell_types"
    } else{ # get scaffold
        scaffold <- spaceRATScaffolds::getScaffold(object)
        return(scaffold)
    }
    # Return subsetted scaffold
    return(list("object" = object, "pheno" = pheno, "colname" = colname))
}
