#' Helper function to load package data
#' @importFrom utils data
#' @noRd
loadData <- function(name){
    utils::data(
        list = name,
        package = "spaceRATScaffolds",
        envir = environment()
        )
    return(get(name))
}

