
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

**Detailed walk-through of how the package works**
