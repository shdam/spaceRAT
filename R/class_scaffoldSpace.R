setOldClass("prcomp")

#' A S4 class to store all the information of scaffold PCA space
#'
#' @slot DEgene a character vector containing names of differentially expressed genes in scaffold dataset
#' @slot label a character vector containing labels of each data point
#' @slot pca a prcomp object returned by function \code{prcomp()}
#' @slot pcs the two indices of principle components to be plot
#' @slot plot_mode a character indicating whether to add tiny label at each data point. Options are "dot" or "tiny_label"
#'
scaffoldSpace <- setClass("scaffoldSpace",
                          slots = list(DEgene="character",
                          label="character",
                          pca="prcomp",
                          pcs="numeric",
                          plot_mode="character")
         )
