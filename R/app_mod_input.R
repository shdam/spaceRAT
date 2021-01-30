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
    # Select space to plot in
    selectInput(
      inputId = ns("space"),
      label = "Select space:",
      choices = c("Scaffold PCA", "gtex", "other"),
      selected = "Scaffold PCA"
    ),
    uiOutput(ns("scaffold")),
    uiOutput(ns("doneScaf")),
    #
    # Input samples
    h2("Select input samples"),
    fileInput(
      inputId = ns('sample_exprs'),
      label = "Choose an expression matrix:",
      multiple = FALSE
    ),
    fileInput(
      inputId = ns('sample_pheno'),
      label = "Choose an phenotype matrix:",
      multiple = FALSE
    ),
    textInput(
      inputId = ns("group"),
      label = "Column name of phenotype data:",
      width = "45%",
      value = "cancer_type"
    ) %>%
      tagAppendAttributes(class = "dim"),
    actionButton(inputId = ns("upSample"),
                 label = "Process"),
    uiOutput(ns("doneSample")),
    # p("Currently, the program only works with the 'processed_data/testdata.Rdata'-file, please upload that file for testing.")%>%
      # tagAppendAttributes(class = "msg"),
    uiOutput(ns("edit")),
    uiOutput(ns("fileCheck")),
    uiOutput(ns("save")),
    uiOutput(ns("plotname"))
  )
}

#' input Server Function
#' @importFrom utils write.csv
#' @importFrom stringr str_c
#' @importFrom readr read_csv
#' @import magrittr dplyr tibble
#' @import shinyWidgets
#' @noRd
mod_input_server <- function(input, output, session, r){
  ns <- session$ns
  # Allow large files
  options(shiny.maxRequestSize = 6000*1024^2) # 6 GB
  # Upload Scaffold files ----
  observeEvent( input$space, {
    if (input$space == "Scaffold PCA"){
      r$column <- "cell_types"
      output$scaffold <- renderUI({
        list(
          # Build Scaffold input
          h2("Build Scaffold"),
          fileInput(
            inputId = ns('scaf_exprs'),
            label = "Choose an expression matrix:",
            multiple = FALSE
          ),
          fileInput(
            inputId = ns('scaf_pheno'),
            label = "Choose an phenotype matrix:",
            multiple = FALSE
          ),
          textInput(
            inputId = ns("column"),
            label = "Column to plot:",
            width = "45%",
            value = r$column
          ) %>%
            tagAppendAttributes(class = "dim"),
          actionButton(inputId = ns("upScaffold"),
                       label = "Process")
        )
        })
    } else {output$scaffold <- NULL}
  })

  # Uploaded scaffold files ----
  observeEvent( input$upScaffold, {
    # Store files
    r$scaf_exprs <- input$scaf_exprs$datapath %>% read.csv()
    r$scaf_pheno <- input$scaf_pheno$datapath %>% read.csv()
    rownames(r$scaf_exprs) <- r$scaf_exprs[,1] # rownames of expression matrix are gene names, colnames are sample names
    r$scaf_exprs <- r$scaf_exprs[,-1]            # every element of matrix should be numbers.
    rownames(r$scaf_pheno) <- r$scaf_pheno[,1]   # rownames of phenotype table are sample names, colnames can be whatever.
    r$scaf_pheno <- r$scaf_pheno[,-1,drop = F]   # drop=F prevents conversion from data.frame to vector
    if(!all(colnames(r$scaf_exprs) == rownames(r$scaf_pheno))) stop("Column names in scaffold expression matrix does not correspond to rows in phonetype table.")
    #
    r$column <- input$column
    output$doneScaf <- renderUI({
      list(
        p("Processed")
        # selectInput(
        #   inputId = ns("source"),
        #   label = "Select data source:",
        #   choices = c("Human", "Mouse"),
        #   selected = "Human"
        # )
      )
    })
    output$edit <- NULL
    output$save <- NULL
  })


  # Uploaded sample files ----
  observeEvent( input$upSample, {
    # Store files
    r$sample_exprs <- input$sample_exprs$datapath %>% read.csv()
    rownames(r$sample_exprs) <- r$sample_exprs[,1] # rownames of expression matrix are gene names, colnames are sample names
    r$sample_exprs <- r$sample_exprs[,-1]            # every element of matrix should be numbers.
    r$sample_pheno <- input$sample_pheno$datapath %>% read.csv()
    rownames(r$sample_pheno) <- r$sample_pheno[,1]   # rownames of phenotype table are sample names, colnames can be whatever.
    r$sample_pheno <- r$sample_pheno[,-1,drop = F]   # drop=F prevents conversion from data.frame to vector
    if(!all(colnames(r$sample_exprs) == rownames(r$sample_pheno))) stop("Column names in sample expression matrix does not correspond to rows in phonetype table.")
    #
    r$group <- input$group
    output$doneSample <- renderUI({
      r$rmPlot <- TRUE
      list(
        p("Processed"),
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
    # Hide stuff
    output$fileCheck <- NULL
    output$doneScaf <- NULL
    output$doneSample <- NULL

    # Input values
    r$validate <- input$validate
    r$source <- input$source
    r$inputFile <- input$inputFile
    r$space <- input$space


    # Plot standard labels
    r$title <- paste("Samples projected on", r$space)
    r$x <- "PC1"
    r$y <- "PC2"
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

