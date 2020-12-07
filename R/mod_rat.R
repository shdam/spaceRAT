#' rat UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom plotly plotlyOutput
mod_rat_ui <- function(id){
  ns <- NS(id)
  tagList(
    plotlyOutput(outputId = ns("pca_plotly"))
  )
}

#' rat Server Function
#' @import magrittr dplyr
#' @importFrom plotly ggplotly renderPlotly layout config add_annotations
#' @noRd
mod_rat_server <- function(input, output, session, r){
  ns <- session$ns
  # Insert load data script ----
  # load(r$inputFile$datapath)
  # mapped <- map_samples(testdata)
  # r$all_classes <- mapped$data$class
  r$classes <- r$all_classes
  # r$data <- mapped$data

  # PCA plot ----
  # pca_plot <- reactive(
  #   mapped$data %>%
  #     filter(class %in% r$classes) %>%
  #     plot_samples() +
  #     labs(title = r$title,
  #          x = r$x,
  #          y = r$y,
  #          subtitle = r$subtitle
  #     )
  # )



  # Render plot ----
  if(FALSE){
  output$pca_plotly <- renderPlotly(
      ggplotly( pca_plot() +
                  theme( legend.title = element_blank() )
              ) %>%
        add_annotations( text="Class", xref="paper", yref="paper",
                         x=1.03, xanchor="left",
                         y=0.85, yanchor="bottom",    # Same y as legend below
                         legendtitle = TRUE, showarrow = FALSE ) %>%
        layout( legend = list(y = 0.8, yanchor = "top" ),
                title = list(text = paste0(r$title,
                                          "<br>",
                                          "<sup>",
                                          r$subtitle,
                                          "</sup>"))) %>%
        config( displayModeBar = FALSE )

  )
  # Save plot ----
  observeEvent(r$savePlot, {
    ggsave(plot = pca_plot(),
           path = getwd(),
           filename = str_c(r$plotname, ".", r$format),
           device = r$format,
           width = r$width,
           height = r$height,
           units = "cm",
           limitsize = TRUE
           )
  })}
}

