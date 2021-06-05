
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RAT

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

## Install as package

### 1\. Create virtual environment with renv

Initialize a local R environment:

``` r
# Open project in Rstudio
# Install and initialize renv 
install.packages("renv")
library(renv)
renv::init()
```

### 2\. Install package from github

``` r
# To ensure Rstudio looks up BioConductor packages run:
setRepositories(ind = c(1:6, 8))
# Then install package with
devtools::install_github("XueningHe/RAT_package")
```

## Install as cloned repository

### 1\. Clone repository

``` sh
# In terminal at desired directory
git clone git@github.com:XueningHe/RAT_package.git
```

### 2\. Install dependecies

``` r
# Open RAT.Rproj in Rstudio

# To ensure Rstudio looks up BioConductor packages run:
setRepositories(ind = c(1:6, 8))
# Install renv
install.packages("renv")
library(renv)
# Restore from lock file
renv::restore()
```

### 3\. Load package

``` r
install.packages("devtools")
library(devtools)
devtools::load_all(getwd())
```

## Usage

### Load Shiny app

``` r
library(RAT)

run_app()
```

*Some text about how the shiny app works.*

### Running without Shiny

**Load external data from directory “/RAT\_package” for illustration.**

The package assume the count matrix to be log transformed. If you want
to input raw data, please specify `data=raw` in function
`buildScaffold()`

``` r
# dmap will be the scaffold dataset

# load expression matrix
dmap_exprs <- read.csv("exprs_dmap.csv") # input data type can be anything like csv  
rownames(dmap_exprs) <- dmap_exprs[,1]   # rownames of expression matrix are gene names, colnames are sample names
dmap_exprs <- dmap_exprs[,-1]            # every element of matrix should be numbers.

# load phenotype table
dmap_pheno <- read.csv("pData_dmap.csv") # input data type can be anything like csv  
rownames(dmap_pheno) <- dmap_pheno[,1]   # rownames of phenotype table are sample names, colnames can be whatever.
dmap_pheno <- dmap_pheno[,-1,drop = F]   # drop=F prevents conversion from data.frame to vector

# check row-col concordance
# colnames of expression matrix should be identical to rownames of phenotype table
all(colnames(dmap_exprs)==rownames(dmap_pheno)) 
#> [1] TRUE

# ilaria will be new samples for projection

# load expression matrix
ilaria_exprs <- read.csv("exprs_ilaria.csv") # input data type can be anything like csv  
rownames(ilaria_exprs) <- ilaria_exprs[,1]   # rownames of expression matrix are gene names, colnames are sample names
ilaria_exprs <- ilaria_exprs[,-1]            # every element of matrix should be numbers.

# load phenotype table
ilaria_pheno <- read.csv("pData_ilaria.csv")  # input data type can be anything like csv  
rownames(ilaria_pheno) <- ilaria_pheno[,1]    # rownames of phenotype table are sample names, colnames can be whatever.
ilaria_pheno <- ilaria_pheno[,-1,drop = F]    # drop=F prevents conversion from data.frame to vector

# check row-col concordance
# colnames of expression matrix should be identical to rownames of phenotype table
all(colnames(ilaria_exprs)==rownames(ilaria_pheno)) 
#> [1] TRUE
```

Users mainly need 2 functions: `buildScaffold` and `projectSample`.

**`buildScaffold()`**

`buildScaffold` takes in an expression matrix (`exprs_scaffold`), a
phenotype table (`pData_scaffold`), a column name within phenotype table
(`group_scaffold`), and returns a scaffoldSpace object, containing a
vector of names of differentially expressed genes in scaffold dataset, a
gg object returned by ggplot, and a prcomp object returned by function .
Type `?buildScaffold` in console to see more optional parameters (which
might be able to be adjusted in the webapp). Type `?scaffoldSpace` in
console to see more information of the S4 class scaffoldSpace defined in
this package.

``` r
library(RAT)
#> Registered S3 methods overwritten by 'ggalt':
#>   method                  from   
#>   grid.draw.absoluteGrob  ggplot2
#>   grobHeight.absoluteGrob ggplot2
#>   grobWidth.absoluteGrob  ggplot2
#>   grobX.absoluteGrob      ggplot2
#>   grobY.absoluteGrob      ggplot2
# To see the scaffold space 
space <- buildScaffold(dmap_exprs,dmap_pheno,"cell_types") # "cell_types" is a column name of dmap_pData
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

For now, calling buildScaffold will automatically print out the scaffold
PCA plot. This can be changed in “R/buildScaffold.R”.

**`projectSample`** `projectSample()` takes in a scaffoldSpace object
output by function `buildScaffold()`, an expression matrix, a
corresponding phenotype table, and a column name. It outputs the
scafffoled PCA plot, together with new samples

``` r
projectSample(space,ilaria_exprs,ilaria_pheno,"cancer_type") # "cancer_type" is a column name of dmap_pData
#> 6 genes are added to count matrix, with imputed expression level 0.
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />
