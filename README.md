
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spaceRAT <img src="inst/spaceRAT.png" width="200" align="right" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

## Installation

*Dependencies*

``` r
# Install Bioconductor dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SummarizedExperiment", "S4Vectors", "BiocStyle"))

# If you run into Matrix/irlba issues, 
# try reinstalling both from source and restart your R session:
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")

# Install spaceRATScaffolds from github
remotes::install_github("shdam/spaceRATScaffolds", build_vignettes = TRUE)
```

*Install from GitHub*

``` r
remotes::install_github("shdam/spaceRAT", build_vignettes = TRUE)
```

## View vignettes

``` r
browseVignettes(package = "spaceRAT")
```

## Usage

It takes two steps to perform ranked analysis of transcriptome:

1.  build a scaffold space
2.  project your new samples onto the scaffold.

Example code will be provided here to illustrate each. Additional
information can be found in the vignettes.

### Build a scaffold space

There are two ways to get a scaffold space. You can either obtain the
prebuilt scaffolds, or build a scaffold space of your own, by providing
a count matrix, a phenotype table, and a column name of the phenotype
table to `buildScaffold()`.

Build a scaffold with example data:

``` r
library(spaceRAT)
#> Loading required package: spaceRATScaffolds
data("DMAP_exprs", "DMAP_pData", package = "spaceRATScaffolds")
scaffold <- buildScaffold(
    DMAP_exprs, pheno = DMAP_pData,
    colname = "cell_types", data = "exprs")
#> Preprocessing complete.
#> Finding differentially expressed genes
#> Reducing dimensions.
#> Scaffold is built.
plotScaffold(
    scaffold, title = "Scaffold Plot",
    dimred = "PCA", dims = c(1,2), plot_mode = "dot")
```

<img src="man/figures/README-build-1.png" width="100%" />

### Project new samples

Get a list of available prebuilt scaffolds with:

``` r
library("spaceRATScaffolds")
listScaffolds()
#> [1] "TCGA.v2" "DMAP.v1" "GTEx.v1" "TCGA.v1"
```

Project a sample of interest into a custom built or prebuilt scaffold:

``` r
# Load count data
data("ilaria_counts", package="spaceRATScaffolds")

# Load custom or prebuilt scaffold
scaffold <- buildScaffold("DMAP") # omitting '.vX' gets the latest version
#> Scaffold will be removed after download.
#> [zen4R][INFO] ZenodoRecord - Download in sequential mode 
#> [zen4R][INFO] ZenodoRecord - Will download 1 file from record '10842509' (doi: '10.5281/zenodo.10842509') - total size: 839 KiB 
#> [zen4R][INFO] Downloading file 'DMAP.v1.rds' - size: 839 KiB
#> [zen4R][INFO] File downloaded at '/net/mimer/mnt/tank/projects2/shd/spaceRAT'.
#> [zen4R][INFO] ZenodoRecord - Verifying file integrity... 
#> [zen4R][INFO] File 'DMAP.v1.rds': integrity verified (md5sum: 1ffc9181d70a414e75458acd4bddb9f2)
#> [zen4R][INFO] ZenodoRecord - End of download
#> Scaffold deleted.

# Project sample
projectSample(
    sample = ilaria_counts,
    scaffold = scaffold, 
    plot_mode = "dot",
    dims = c(1,2),
    title = "Samples projected into DMAP scaffold in PCA space")
#> Preprocessing complete.
#> 6 genes are added to count matrix with imputed expression level 0.
```

<img src="man/figures/README-project-1.png" width="100%" />
