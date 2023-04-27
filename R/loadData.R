#' Helper function to load package data
#' @noRd
loadData <- function(name){
    data_file <- system.file("data", paste0(name,".rda"), package = "spaceRAT")
    load(data_file)
    return(get(name))
}
