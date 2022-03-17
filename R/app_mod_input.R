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
      choices = c("dmap", "gtex", "other"), #"Build scaffold"
      selected = "dmap"
    ),
    uiOutput(ns("scaffold")),
    uiOutput(ns("doneScaf")),
    #
    # Input samples
    # h2("Select input samples"),
    uiOutput( # Expression tooltip
      outputId = ns("exprinfo")
    ),
    fileInput(
      inputId = ns('sample_exprs'),
      label = "Choose an expression matrix:",
      multiple = FALSE
    ),
    uiOutput( # Phenotype tooltip
      outputId = ns("phenoinfo")
    ),
    fileInput(
      inputId = ns('sample_pheno'),
      label = "Choose a phenotype matrix:",
      multiple = FALSE
    ),
    checkboxInput(
      inputId = ns("loading"),
      label = "Loading plot",
      value = FALSE
    ),
    selectInput(
      inputId = ns("plot_mode"),
      label = "Plot mode",
      choices = c("dot", "tiny_label"),
      selected = "dot",
      width = "45%"
    ),
    uiOutput(ns("group")),
    # actionButton(inputId = ns("upSample"),
    #              label = "Process"),
    uiOutput(ns("classes")),
    actionButton(
        inputId = ns("validate"),
        label = "Plot"
    ),
    uiOutput(ns("noFiles")),
    # p("Currently, the program only works with the 'processed_data/testdata.Rdata'-file, please upload that file for testing.")%>%
      # tagAppendAttributes(class = "msg"),
    uiOutput(ns("edit")),
    uiOutput(ns("fileCheck")),
    br(),
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
  options(shiny.maxRequestSize = 6000*1000^2) # 6 GB
  # Upload Scaffold files ----
  observeEvent( input$space, {
    r$space <- input$space
    # Add inputs related to building your own scaffold
    if (input$space == "Build scaffold"){
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

    # Initiate values ----
    if(input$space == "dmap"){
      output$noFiles <- NULL
      r$all_classes <- unique(pData_dmap$cell_types)
      r$classes <- r$all_classes
      r$scaf_exprs <- exprs_dmap
      r$scaf_pheno <- pData_dmap
    } else if(input$space == "gtex"){
      # r$inputFile <- r$
      output$noFiles <- renderUI({
        p("Not yet implemented")
      })
    }else if(input$space == "other"){
      # r$inputFile <- r$
      output$noFiles <- renderUI({
        p("Not yet implemented")
      })}
      #
    output$classes <- renderUI(list(
      selectInput(
        inputId = ns("classes"),
        label = "Select cell types to include",
        choices = sort(r$all_classes),
        multiple = TRUE,
        selected = r$classes
        ),
      p("Mark and press the 'delete' button to remove classes.") %>%
        tagAppendAttributes(class = 'msg')
    ))
  })


  observeEvent( r$sample_pheno, {
    #   r$plot_mode <- input$plot_mode
    #
    if(!(r$column %in% colnames(r$sample_pheno))) r$column <- colnames(r$sample_pheno)
    if(!is.null(r$sample_pheno)){
      output$group <- renderUI({
        selectInput(
          inputId = ns("group"),
          label = "Column name in phenotype:",
          width = "45%",
          choices = colnames(r$sample_pheno),
          selected = colnames(r$sample_pheno)[-1]
        ) %>%
          tagAppendAttributes(class = "dim")
      })
    } else{
      output$group <- NULL
    }

  })


  # Tooltips ----
  tooltip_expr <- tibble::tribble(
    ~` `, ~sample1, ~sample2, ~`...`,
    "gene1", 200, 900, "...",
    "gene2", 1300, 5, "..."
  )
  tooltip_pheno <- tibble::tribble(
    ~` `, ~cell_types,
    "sample1", "type1",
    "sample2", "type2",
    "...", "..."
  )
  output$exprinfo <- renderUI({
      shinyWidgets::dropdownButton(
        h4("Example expression matrix layout"),
        renderTable(tooltip_expr),
        inline = TRUE,
        size = "xs",
        circle = TRUE,
        # status = "danger",
        icon = icon("info"), width = "500px",
        tooltip = shinyWidgets::tooltipOptions(title = "Click to see expression matrix example")
      )
  })
  output$phenoinfo <- renderUI({
      shinyWidgets::dropdownButton(
        h4("Example phenotype matrix layout"),
        renderTable(tooltip_pheno),
        inline = TRUE,
        size = "xs",
        circle = TRUE,
        # status = "danger",
        icon = icon("info"), width = "500px",
        tooltip = shinyWidgets::tooltipOptions(title = "Click to see phenotype matrix example")
      )
  })

  # Uploaded sample files ----
  observeEvent(input$sample_exprs, {
    output$noFiles <- NULL
    output$edit <- NULL
    output$save <- NULL
    r$rmPlot <- TRUE
    if(!is.null(input$sample_exprs)) r$exprs_datapath <- input$sample_exprs$datapath
  })
  observeEvent(input$sample_pheno, {
    if(!is.null(input$sample_pheno)) r$pheno_datapath <- input$sample_pheno$datapath
  })

  # Plot button is pressed ----
  observeEvent( input$validate, {
    r$loading <- input$loading
    # Hide stuff
    output$fileCheck <- NULL
    output$doneScaf <- NULL

    if(input$space == "Build scaffold"){
      # Store scaffold files
      r$scaf_exprs <- input$scaf_exprs$datapath %>% read.csv()
      r$scaf_pheno <- input$scaf_pheno$datapath %>% read.csv()

    }
    # If no input
    if(input$space != "Build scaffold" & (is.null(input$scaf_exprs) & is.null(input$scaf_pheno))){
    if(!is.null(input$sample_exprs)){ # & !is.null(input$sample_pheno)

    if(is.null(input$sample_pheno)) r$sample_pheno <- NULL

    # Input values
    r$column <- input$group
    r$plot_mode <- input$plot_mode
    r$space <- input$space
    r$validate <- input$validate
    r$classes <- input$classes


    # Dropdown edit button ----
    # Sys.sleep(3)
    output$edit <- renderUI({
      list(
        dropdown(
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
          # textInput(
          #   inputId = ns("x"),
          #   label = "x",
          #   width = "45%",
          #   value = r$x
          # ) %>%
          #   tagAppendAttributes(class = "dim"),
          # textInput(
          #   inputId = ns("y"),
          #   label = "y",
          #   width = "45%",
          #   value = r$y
          # ) %>%
          #   tagAppendAttributes(class = "dim"),

          # Dropdown styling
          style = "unite", icon = icon("gear"),
          color = "danger",
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
        dropMenu(
          actionButton(ns("down"), "Download"),
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
      ))
    })
    } else{ # End if no input ----
      output$noFiles <- renderUI(
        p("A file is missing!")
      )
    }
    } else{ # No scaffold given
      output$noFiles <- renderUI(
        p("A file is missing!")
      )
    }
  })

  # Apply changes ----
  observeEvent( input$applyBtn, {
    r$title <- input$title
    r$subtitle <- input$subtitle
    output$plotname <- NULL
    r$plotname <- NULL
    # r$validate <- input$applyBtn
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

