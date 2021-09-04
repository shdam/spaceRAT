#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny shinythemes
#' @noRd
app_ui <- function(request) {
        tagList(
                # Leave this function for adding external resources
                golem_add_external_resources(),
                # List the first level UI elements here
                fluidPage(
                        # theme = shinytheme("spacelab"),
                        headerPanel(
                                div(class = "header",
                                    h1(class = "header", "RAT")
                                    )
                                ),
                        sidebarLayout(
                                sidebarPanel(
                                        mod_input_ui("input_ui_1")
                                ),
                                mainPanel(
                                        mod_rat_ui("rat_ui_1")
                                )
                        )


                )
        )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){

        add_resource_path(
                'www', app_sys('app/www')
        )

        tags$head(
                favicon(ext = "png"),
                bundle_resources(
                        path = app_sys('app/www'),
                        app_title = 'RAT'
                ),
                # Add here other external resources
                # for example, you can add shinyalert::useShinyalert()
                tags$link(rel="stylesheet", type="text/css", href="www/custom.css")
        )
}

