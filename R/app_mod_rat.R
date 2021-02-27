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
#' @importFrom plotly ggplotly renderPlotly layout config add_annotations
#' @noRd
mod_rat_server <- function(input, output, session, r){
  ns <- session$ns

  # Prepare scaffold
  # rownames(r$scaf_exprs) <- r$scaf_exprs[,1] # rownames of expression matrix are gene names, colnames are sample names
  # r$scaf_exprs <- r$scaf_exprs[,-1]            # every element of matrix should be numbers.
  # rownames(r$scaf_pheno) <- r$scaf_pheno[,1]   # rownames of phenotype table are sample names, colnames can be whatever.
  # r$scaf_pheno <- r$scaf_pheno[,-1,drop = F]   # drop=F prevents conversion from data.frame to vector
  # if(!all(colnames(r$scaf_exprs) == rownames(r$scaf_pheno))) stop("Column names in scaffold expression matrix does not correspond to rows in phonetype table.")
  #



  # Prepare sample
  # Store sample files
  r$sample_exprs <- r$exprs_datapath %>% read.csv()
  rownames(r$sample_exprs) <- r$sample_exprs[,1] # rownames of expression matrix are gene names, colnames are sample names
  r$sample_exprs <- r$sample_exprs[,-1]            # every element of matrix should be numbers.
  r$sample_pheno <- r$pheno_datapath %>% read.csv()
  rownames(r$sample_pheno) <- r$sample_pheno[,1]   # rownames of phenotype table are sample names, colnames can be whatever.
  r$sample_pheno <- r$sample_pheno[,-1,drop = F]   # drop=F prevents conversion from data.frame to vector
  if(!all(colnames(r$sample_exprs) == rownames(r$sample_pheno))) stop("Column names in sample expression matrix does not correspond to rows in phonetype table.")
  #

  r$all_classes <- unique(r$scaf_pheno[,r$column])
  r$classes <- r$all_classes

  r$sample_exprs_fin <- r$sample_exprs
  r$sample_pheno_fin <- r$sample_pheno
  # group <- r$group



  # Identify space to use
  g <- reactive(
    buildScaffold(r$scaf_exprs, r$scaf_pheno, r$column, r$classes)
  )

  pca_plot <- reactive(
    projectSample(space = g(),
                  exprs_sample = r$sample_exprs_fin,
                  pData_sample = r$sample_pheno_fin,
                  group_sample = r$group,
                  title = r$title,
                  y = r$y,
                  x = r$x)
  )


  # Insert load data script ----
  # load(r$inputFile$datapath)
  # mapped <- map_samples(testdata)
  # r$all_classes <- mapped$data$class
  # r$classes <- r$all_classes
  # r$data <- mapped$data

  # PCA plot ----
  # pca_plot <- reactive(
    # mapped$data %>%
    #   filter(class %in% r$classes) %>%
    #   plot_samples() +
    #   labs(title = r$title,
    #        x = r$x,
    #        y = r$y,
    #        subtitle = r$subtitle
    #   )
  # )

  # eset_dmap <- createEset(expr_mat = exprs_dmap,
  #                         pheno = pData_dmap)
  #
  # DEgenes <- findDEGenes(eset = eset_dmap,
  #                        group = "cell_types",
  #                        pval_cutoff = 0.05,
  #                        lfc_cutoff = 2)
  #
  # g <- buildScaffold(exprs_scaffold = exprs_dmap,
  #                    pData_scaffold = pData_dmap,
  #                    group_scaffold = "cell_types")
  #
  # pca_plot <- reactive(
  #   projectSample(space = g,
  #                 exprs_sample = exprs_ilaria,
  #                 pData_sample = pData_ilaria,
  #                 group_sample = "cancer_type",
  #                 title = r$title)
  # )

  # Render plot ----
  output$pca_plotly <- renderPlotly(

    ggplotly(
      pca_plot() +
        theme( legend.title = element_blank())
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
      config(displayModeBar = FALSE)
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
           limitsize = FALSE
           )
  })
}

