#' Helper function to load package data
#' @noRd
loadData <- function(name){
    data(list = name, package = "spaceRATScaffolds")
    return(get(name))
}
