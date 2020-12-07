#' input UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_input_ui <- function(id){
  ns <- NS(id)
  tagList(
    fileInput(
      inputId = ns('inputFile'),
      label = "Choose an input file:",
      multiple = TRUE
    ),
    p("Currently, the program only works with the 'processed_data/testdata.Rdata'-file, please upload that file for testing.")%>%
      tagAppendAttributes(class = "msg"),
    uiOutput(ns("edit")),
    uiOutput(ns("fileCheck")),
    uiOutput(ns("save")),
    uiOutput(ns("plotname"))
  )
}

#' input Server Function
#' @importFrom utils write.csv
#' @importFrom stringr str_c
#' @import shinyWidgets
#' @noRd
mod_input_server <- function(input, output, session, r){
  ns <- session$ns
  # Allow large files
  options(shiny.maxRequestSize = 6000*1024^2) # 6 GB
  # Uploaded files ----
  observeEvent( input$inputFile, {
    output$fileCheck <- renderUI({
      list(
        selectInput(
          inputId = ns("source"),
          label = "Select data source:",
          choices = c("Human", "Mouse"),
          selected = "Human"
        ),
        selectInput(
          inputId = ns("space"),
          label = "Select space:",
          choices = c("gtex", "other"),
          selected = "gtex"
        ),
        actionButton(
          inputId = ns("validate"),
          label = "Plot"
        )
      )
    })
    output$edit <- NULL
    output$save <- NULL
  })
  # Plot button is pressed ----
  observeEvent( input$validate, {
    # Remove option to change source
    output$fileCheck <- NULL

    # Input values
    r$validate <- input$validate
    r$source <- input$source
    r$inputFile <- input$inputFile
    r$space <- input$space

    # Plot standard labels
    r$title <- paste("My data in", r$space, "space")
    r$x <- "PC5"
    r$y <- "PC9"
    r$subtitle <- str_c("Data source: ", input$source)

    # Dropdown edit button ----
    output$edit <- renderUI({
      list(
        dropdown(
          selectInput(
            inputId = ns("classes"),
            label = "Select classes to plot",
            choices = sort(r$all_classes),
            multiple = TRUE,
            selected = r$classes
          ),
          p("Mark and press the 'delete' button to remove classes.") %>%
            tagAppendAttributes(class = 'msg'),
          textInput(
            inputId = ns("title"),
            label = "Title",
            value = r$title
          ),
          textInput(
            inputId = ns("subtitle"),
            label = "Subtitle",
            value = r$subtitle
          ),
          textInput(
            inputId = ns("x"),
            label = "x",
            width = "45%",
            value = r$x
          ) %>%
            tagAppendAttributes(class = "dim"),
          textInput(
            inputId = ns("y"),
            label = "y",
            width = "45%",
            value = r$y
          ) %>%
            tagAppendAttributes(class = "dim"),

          # Dropdown styling
          style = "unite", icon = icon("gear"),
          status = "primary", width = "300px",
          tooltip = tooltipOptions(title = "Click to edit plot"),
          animate = animateOptions(
            enter = shinyWidgets::animations$fading_entrances$fadeInLeftBig,
            exit = shinyWidgets::animations$fading_exits$fadeOutRightBig
          ),
          actionButton(
            inputId = ns("applyBtn"),
            label = "Apply changes"
          )
        )
      )
    })
    output$plotname <- NULL
    r$plotname <- NULL

    # Save plot ----
    output$save <- renderUI({
      list(
        radioButtons(
          inputId = ns("format"),
          label = "Choose output format:",
          choices = c("png", "jpeg", "pdf", "svg"),
          selected = "png",
          inline = TRUE
        ),
        numericInput(
          inputId = ns("width"),
          label = "Plot width (cm)",
          value = 20,
          min = 5,
          width = "45%",
          step = 0.5
        ) %>%
          tagAppendAttributes(class = "dim"),
        numericInput(
          inputId = ns("height"),
          label = "Plot height (cm)",
          value = 12,
          min = 5,
          width = "45%",
          step = 0.5
        ) %>%
          tagAppendAttributes(class = "dim"),
        actionButton(
          inputId = ns("savePlot"),
          label = "Download plot"
        ),
        # Download button ----
        downloadButton(ns("downloadData"), "Download data")
      )
    })
  })

  # Apply changes ----
  observeEvent( input$applyBtn, {
    r$title <- input$title
    r$x <- input$x
    r$y <- input$y
    r$classes <- input$classes
    r$source <- input$source
    r$subtitle <- input$subtitle
    output$plotname <- NULL
    r$plotname <- NULL
    })
  # Press save ----
  observeEvent( input$savePlot, {
    r$savePlot <- input$savePlot
    r$format <- input$format
    r$height <- input$height
    r$width <- input$width
    r$plotname <- "plot_rat"
  })
  # Print filename ----
  observeEvent( r$plotname, {
    if ( !is.null(r$plotname) ){
      output$plotname <- renderUI(
        p(str_c("Plot saved with filename '", str_c(r$plotname, ".", r$format), "'"))
      )
    }
  })
  # Download data ----
  output$downloadData <- downloadHandler(
    filename = "RAT_data.csv",
    content = function(file) {
      write.csv(r$data, file, row.names = FALSE)
    }
  )
}

