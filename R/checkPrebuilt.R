#' Check if input is a prebuilt scaffold
#'
#' @noRd
checkPrebuilt <- function(object, classes){
    if(
        is(object, "character") &&
        object %in% listScaffolds() &&
        is(classes, "NULL")
    ){ # Return prebuil scaffold
        space <- loadData(object)
        return(space)
        # prebuilt DMAP samples removed
    } else if(
        is(object, "character") &&
        object == "DMAP_scaffold" &&
        !is(classes, "NULL")
    ){ # Subsetting a Scaffold
        object <- loadData("exprs_dmap")
        pheno <- loadData("pData_dmap")
        colname <- "cell_types"
    } else if(is(object, "character")){ # Unknown scaffold
        stop(
            "Incorrectly named prebuilt scaffold. The available are: ",
            listScaffolds()
        )
    }
    # Return subsetted scaffold
    return(list("object" = object, "pheno" = pheno, "colname" = colname))
}
