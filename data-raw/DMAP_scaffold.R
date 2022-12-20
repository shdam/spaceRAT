## code to prepare `DMAP_scaffold` dataset goes here

library(spaceRAT)

counts <- readr::read_csv("tmp/exprs_dmap.csv", show_col_types = FALSE)
pdata <- readr::read_csv("tmp/pData_dmap.csv", show_col_types = FALSE)
DMAP_scaffold <- buildScaffold(
    object = counts,
    pheno_scaffold = pdata,
    colname = "cell_types"
)


usethis::use_data(DMAP_scaffold, overwrite = TRUE)
