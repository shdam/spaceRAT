#' A S4 class to store all the information of scaffold PCA space
#' @slot DEgenes a vector containing names of differentially expressed genes in scaffold dataset
#' @slot graph a gg object returned by ggplot
#' @slot pca a prcomp object returned by function \code{prcomp()}
#'
setOldClass("gg")
setOldClass("prcomp")
scaffoldSpace <- setClass("scaffoldSpace",
                          slots = list(DEgenes="character",
                          graph="gg",
                          pca="prcomp")
         )
