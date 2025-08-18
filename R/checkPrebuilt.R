#' Check if input is a prebuilt scaffold
#'
#' @noRd
checkPrebuilt <- function(object, classes, path){
    stopifnot("Please only check a character string" = is(object, "character"))
    if( # Return prebuilt scaffold
        object == "DMAP.v1" & !is.null(classes)
        ){
        object <- loadData("DMAP_exprs")
        pheno <- loadData("DMAP_pData")
        colname <- "cell_types"
    } else{ # get scaffold
        scaffold <- spaceRATScaffolds::getScaffold(
          object, store = !is.null(path), path = ifelse(is.null(path), "scaffolds", "path"))
        return(scaffold)
    }
    # Return subsetted scaffold
    return(list("object" = object, "pheno" = pheno, "colname" = colname))
}
