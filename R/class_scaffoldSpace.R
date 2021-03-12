setOldClass("gg")
setOldClass("prcomp")

#' A S4 class to store all the information of scaffold PCA space
#'
#' @slot DEgenes a vector containing names of differentially expressed genes in scaffold dataset
#' @slot graph a gg object returned by ggplot
#' @slot pca a prcomp object returned by function \code{prcomp()}
#' @slot pcs the two indices of principle components to be plot, by default pcs=c(1,2)
#'
scaffoldSpace <- setClass("scaffoldSpace",
                          slots = list(DEgenes="character",
                          graph="gg",
                          pca="prcomp",
                          pcs="numeric")
         )
