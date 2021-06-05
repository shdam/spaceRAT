
#' Creates PCA plot of test data on the gtex scaffold
#' @param input_data The data to be mapped onto defined scaffold space
#' @import dplyr magrittr tibble
#' @importFrom mclust adjustedRandIndex
#' @importFrom stats kmeans prcomp
#' @noRd
map_samples <- function(input_data){
  ## Load scaffold
  # load(file="./R/sysdata.rda")
  ## Import expression from samples (this will be in any type of format - pre-processing script will follow)
  # testdata <- input_data %>%
    # as_tibble(rownames = "genes") #load("./processed_data/testdata.Rdata")

  ## Subset both sets based on genes in common
  gtex_testdata <- gtex_testdata %>%
    semi_join(gtex_scaffold, by = "genes") %>%
    arrange(.data$genes)

  gtex_scaffold <- gtex_scaffold %>%
    semi_join(gtex_testdata, by = "genes") %>%
    arrange(.data$genes)

  ## Rank both

  gtex_scaffold_rank <- gtex_scaffold %>%
    column_to_rownames("genes") %>%
    mutate_all(.funs = rank)

  testdata_rank <- gtex_testdata %>%
    column_to_rownames("genes") %>%
    mutate_all(.funs = rank)

  ## PCA of Xena
  pca <- gtex_scaffold_rank %>%
    t() %>%
    prcomp(scale. = TRUE)


  ## Project into PC space
  pca_new <- testdata_rank %>%
    t() %>%
    scale(center = pca$center,
          scale = pca$scale) %*%
    pca$rotation

  # Pick PCs based on k means and adjusted rand index
  tests <- expand.grid(c(1:10), c(1:10)) %>%
    filter(.data$Var2 > .data$Var1)
  ARI <- c()
  for(i in 1:nrow(tests)) {
    set.seed(42)
    predicted <- kmeans(x = cbind(pca$x[, tests[i, 1]],
                                  pca$x[, tests[i, 2]]),
                        centers = length(unique(class_sample)),
                        nstart = 25)
    ARI <- append(ARI, adjustedRandIndex(class_sample,
                                         unname(predicted$cluster)))
  }
  d1 <- tests[which.max(ARI), 1]
  d2 <- tests[which.max(ARI), 2]

  ## Plot (should be updated for niceness - something like repelled labels to centroids instead of legend)
  data <- tibble(dim1 = append(pca$x[, d1], pca_new[, d1]),
                 dim2 = append(pca$x[, d2], pca_new[, d2]),
                 class = append(class_sample, rep("newdata", ncol(gtex_testdata %>% select(-.data$genes)))))


  return(list("data" = data, "d1" = d1, "d2" = d2))
}

#' Function that makes a ggplot from the output of "map_samples"
#' @param data data object from the output of "map_samples" containing the two principal components from the analysis
#' @param d1 principal component number of the x-axis
#' @param d2 principal component number of the y-axis
#' @import ggplot2
#' @importFrom dplyr case_when
#' @noRd
plot_samples <- function(data, d1 = 5, d2 = 9){
  pca_plot <- data %>%
    mutate(shape = case_when(class == "newdata" ~ 4,
           TRUE ~ 1) %>% as.factor()) %>%
    ggplot(aes(x = .data$dim1,
               y = .data$dim2,
               color = .data$class)) +
      geom_point(aes(shape = shape)) +
      ggplot2::scale_shape_manual(values = c(1, 4)) +
      labs(x = str_c("PC", d1),
           y = str_c("PC", d2),
           color = "Class"
          )
  return(pca_plot)
}



