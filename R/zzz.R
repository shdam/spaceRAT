#' This function unloads the spaceRATScaffolds package namespace
#' when unloading spaceRAT package
#' @noRd
.onUnload <- function(...) {
    ns <- getNamespace("spaceRATScaffolds")
    if (!is.null(ns)) {
        try(unloadNamespace("spaceRATScaffolds"), silent = TRUE)
    }
}
