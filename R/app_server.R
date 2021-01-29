#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @param r a `reactiveValues()`list with a `session` and `source` element in it.
#' The `r$session` is used to refer to the main session for tab manipulation within modules
#' `r$source` specifies whether the data is from human or a mouse
#' @import shiny
#' @noRd
app_server <- function( input, output, session ) {
        r <- reactiveValues(warning = FALSE)
        r$session <- session

        # Stop app when window is closed
        session$onSessionEnded(function() {
                stopApp()
        })

        # List the first level callModules here
        callModule(mod_input_server, "input_ui_1", r)
        observeEvent(r$validate, {
                callModule(mod_rat_server, "rat_ui_1", r)
        })
}
