## code to prepare `DMAP_scaffold` dataset goes here

library(RAT)

counts <- readr::read_csv("data-raw/exprs_dmap.csv")
pdata <- readr::read_csv("data-raw/pData_dmap.csv")
scaffold <- buildScaffold(
  counts_scaffold = counts,
  pheno_scaffold = pdata,
  colname = "cell_types"
)


usethis::use_data(DMAP_scaffold, overwrite = TRUE)
