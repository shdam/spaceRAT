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
    uiOutput(ns("error")),
    plotOutput(outputId = ns("pca_plot"))
  )
}

#' rat Server Function
#' @importFrom plotly ggplotly renderPlotly layout config add_annotations
#' @importFrom ggalt geom_encircle
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


## Initial plot
  observeEvent( r$space, {
    # Loading plot
    if(r$space == "dmap"){
      output$pca_plot <- renderPlot(height = 700, {
        buildScaffold("prebuilt_DMAP",
                      auto_plot = FALSE,
                      classes = NULL,
                      plot_mode = r$plot_mode) %>%
          plotScaffold(title = r$title)
      })




    } else if(r$space == "gtex"){
      output$pca_plot <- NULL
    } else{
      output$pca_plot <- NULL
    }
  })

  observeEvent( r$exprs_datapath, {
    # Prepare sample
    # Store sample files
    r$sample_exprs <- r$exprs_datapath %>% read.csv(row.names = 1)
    # rownames(sample_exprs) <- sample_exprs[,1] # rownames of expression matrix are gene names, colnames are sample names
    # r$sample_exprs <- sample_exprs[,-1]            # every element of matrix should be numbers.
    # rm(sample_exprs)
  })
  observeEvent( r$pheno_datapath, {
    r$sample_pheno <- r$pheno_datapath %>% read.csv(row.names = 1)
    # rownames(sample_pheno) <- sample_pheno[,1]   # rownames of phenotype table are sample names, colnames can be whatever.
    # r$sample_pheno <- sample_pheno[,-1,drop = F]   # drop=F prevents conversion from data.frame to vector
    # rm(sample_pheno)
  })

  observeEvent( r$validate, {
    if(r$space == "dmap"){
      if(!is.null(r$pheno_datapath)){
        if(!all(colnames(r$sample_exprs) == rownames(r$sample_pheno))) {
          output$error <- renderUI(p("Column names in sample expression matrix does not correspond to rows in the phenotype table.", id = "error"))
          output$pca_plot <- NULL
          stp <- TRUE
        } else if(!(r$column %in% colnames(r$sample_pheno))){
          output$error <- renderUI(p(paste0("A column of the name: '", r$column,"' was not found in the provided phenotype data. Please check that the column names match."), id = "error"))
          output$pca_plot <- NULL
          stp <- TRUE
        } else{
          stp <- FALSE
        }
      } else{
        r$sample_pheno <- NULL
        stp <- FALSE
      }

      #
      if(!stp){


        # Create scaffold
        if(all(r$all_classes %in% r$classes)) classes <- NULL else classes <- r$classes
        r$g <- buildScaffold("prebuilt_DMAP",
                        auto_plot = FALSE,
                        classes = classes,
                        plot_mode = r$plot_mode)





        # r$all_classes <- unique(r$sample_pheno[,r$column])
        # r$classes <- r$all_classes
        # group <- r$group

        # Identify space to use
        # g <- reactive(
        #   buildScaffold(r$scaf_exprs, r$scaf_pheno, r$column, r$classes)
        # )

        pca_plot <- projectSample(space = r$g,
                        counts_sample = r$sample_exprs,
                        pheno_sample = r$sample_pheno,
                        colname = r$column,
                        # classes = r$classes,
                        # group_sample = r$group,
                        title = r$title
                        )
          # Add loading plot if needed
          if(r$loading){
            df <- loadingPlot(r$g, df_only = TRUE)
            pca_plot <- pca_plot +
              ggplot2::geom_text(data = df, ggplot2::aes_string(colnames(df)[2], colnames(df)[3],label=colnames(df)[1],color=colnames(df)[4]),size = 3, angle=30,show.legend = FALSE)+
              ggplot2::geom_segment(data = df, ggplot2::aes_string(x=0, y=0, xend=colnames(df)[2], yend=colnames(df)[3],color=colnames(df)[4]), arrow=ggplot2::arrow(length=ggplot2::unit(0.1,"cm")), alpha=0.75) +
              ggplot2::xlim(min(df[2])*1.5,max(df[2])*1.5)+
              ggplot2::ylim(min(df[3])*1.5,max(df[3])*1.5)
          }

      }

    } else if(r$space == "gtex"){
      # load(r$inputFile$datapath)
      mapped <- map_samples(gtex_testdata)
      # r$all_classes <- mapped$data$class
      # r$classes <- r$all_classes
      # r$data <- mapped$data

      # PCA plot ----
      pca_plot <- reactive(
        mapped$data %>%
          # dplyr::filter(class %in% r$classes) %>%
          plot_samples() +
          labs(title = r$title,
               x = "PC1",
               y = "PC2",
               subtitle = r$subtitle
          )
      )
    }



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
    pca_plot_reactive <- reactive(
      pca_plot +
        # geom_encircle(aes(group=class, fill=class),alpha=0.4) +
        ggplot2::theme_bw() +
        labs(title = r$title,
             subtitle = r$subtitle) +
        theme( legend.title = element_blank())
    )
  output$pca_plot <- renderPlot(height = 700,
                                {pca_plot_reactive()})

  })
  # Save plot ----
  observeEvent(r$savePlot, {
    ggsave(plot = pca_plot_reactive(),
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

