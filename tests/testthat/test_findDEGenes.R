data("DMAP_exprs", package = "spaceRATScaffolds")
data("DMAP_pData", package = "spaceRATScaffolds")
object <- preprocess(
    DMAP_exprs,
    data = "exprs",
    colname = "cell_types",
    pheno = DMAP_pData)

counts_scaffold <- object[[1]]
cell_types <- object[[2]]
rm(object)

test_that("findDEGenes() returns a vector of gene names",{

    DEgenes <- findDEGenes(counts_scaffold, cell_types$cell_types, 0.05,2)
    expect_type(DEgenes,"character")
    expect_true(all(DEgenes %in% rownames(counts_scaffold)))
})

# Test that the function returns an empty vector when no genes pass the thresholds
test_that("findDEGenes returns empty vector when no genes pass thresholds", {
    DEgenes <- spaceRAT:::findDEGenes(counts_scaffold, cell_types$cell_types, pval_cutoff = 1, lfc_cutoff = 100)
    expect_equal(length(DEgenes), 0)
})
