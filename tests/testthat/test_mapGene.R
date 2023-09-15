data("gene_id_converter_hs", package = "spaceRATScaffolds")
test_that("mapGene() converts from Ensembl_gene to others",{
    ensembl <- rownames(exprs_dmap)
    df <- mapGene(ensembl,to="ensembl_transcript")
    tx <- as.character(df[,2])
    expect_match(tx,"^^ENST[0-9]+$")
    expect_false(any(is.na(tx)))
    expect_false("" %in% tx)

    df <- mapGene(ensembl,to="entrez")
    entrez <- as.character(df[,2])
    expect_match(entrez,"^[0-9]+$")
    expect_false(any(is.na(entrez)))
    expect_false("" %in% entrez)

    df <- mapGene(ensembl,to="hgnc_symbol")
    hgnc <- as.character(df[,2])
    expect_true(all(hgnc %in% gene_id_converter_hs$hgnc_symbol))
    expect_false(any(is.na(hgnc)))
    expect_false("" %in% hgnc)

    df <- mapGene(ensembl,to="refseq_mrna")
    refseq <- as.character(df[,2])
    expect_match(refseq,"^NM_")
    expect_false(any(is.na(refseq)))
    expect_false("" %in% refseq)
})


test_that("mapGene() converts from Ensembl_transcript to others",{
    ensembl <- rownames(exprs_dmap)
    df <- mapGene(ensembl,to="ensembl_transcript")
    tx<- df[,2]

    df <- mapGene(tx,to="ensembl_gene")
    ensembl <- df[,2]
    expect_match(ensembl,"^ENSG[0-9]+$")
    expect_false(any(is.na(ensembl)))
    expect_false("" %in% ensembl)

    df <- mapGene(tx,to="entrez")
    entrez <- as.character(df[,2])
    expect_match(entrez,"^[0-9]+$")
    expect_false(any(is.na(entrez)))
    expect_false("" %in% entrez)

    df <- mapGene(tx,to="hgnc_symbol")
    hgnc <- as.character(df[,2])
    expect_true(all(hgnc %in% gene_id_converter_hs$hgnc_symbol))
    expect_false(any(is.na(hgnc)))
    expect_false("" %in% hgnc)

    df <- mapGene(tx,to="refseq_mrna")
    refseq <- as.character(df[,2])
    expect_match(refseq,"^NM_")
    expect_false(any(is.na(refseq)))
    expect_false("" %in% refseq)
})


test_that("mapGene() converts from Entrez to others",{
  ensembl <- rownames(exprs_dmap)[1:5]
  df <- mapGene(ensembl,to="entrez")
  entrez <- as.character(df[,2])

  df <- mapGene(entrez,to="ensembl_gene")
  ensembl <- df[,2]
  expect_match(ensembl,"^ENSG[0-9]+$")
  expect_false(any(is.na(ensembl)))
  expect_false("" %in% ensembl)

  df <- mapGene(entrez,to="ensembl_transcript")
  tx <- df[,2]
  expect_match(tx,"^ENST[0-9]+$")
  expect_false(any(is.na(tx)))
  expect_false("" %in% tx)

  df <- mapGene(entrez,to="hgnc_symbol")
  hgnc <- df[,2]
  expect_equal(length(hgnc),5)
  expect_true(all(hgnc %in% gene_id_converter_hs$hgnc_symbol))
  expect_false(any(is.na(hgnc)))
  expect_false("" %in% hgnc)

  df <- mapGene(entrez,to="refseq_mrna")
  refseq <- df[,2]
  expect_match(refseq,"^NM_")
  expect_false(any(is.na(refseq)))
  expect_false("" %in% refseq)
})


test_that("mapGene() converts between hgnc_symbol to others",{
  ensembl <- rownames(exprs_dmap)
  df <- mapGene(ensembl,to="hgnc_symbol")
  hgnc <- df[,2]

  df <- mapGene(hgnc,to="ensembl_gene")
  ensembl <- df[,2]
  expect_match(ensembl,"^ENSG[0-9]+$")
  expect_false(any(is.na(ensembl)))
  expect_false("" %in% ensembl)

  df <- mapGene(hgnc,to="ensembl_transcript")
  tx <- df[,2]
  expect_match(tx,"^ENST[0-9]+$")
  expect_false(any(is.na(tx)))
  expect_false("" %in% tx)

  df <- mapGene(hgnc,to="entrez")
  entrez <- as.character(df[,2])
  expect_match(entrez,"^[0-9]+$")
  expect_false(any(is.na(entrez)))
  expect_false("" %in% entrez)

  df <- mapGene(hgnc,to="refseq_mrna")
  refseq <- df[,2]
  expect_match(refseq,"^NM_")
  expect_false(any(is.na(refseq)))
  expect_false("" %in% refseq)
})


test_that("mapGene() converts from refseq_mrna to others",{
  ensembl <- rownames(exprs_dmap)
  df <- mapGene(ensembl,to="refseq_mrna")
  refseq <- df[,2]

  df <- mapGene(refseq,to="ensembl_gene")
  ensembl <- df[,2]
  expect_match(ensembl,"^ENSG[0-9]+$")
  expect_false(any(is.na(ensembl)))
  expect_false("" %in% ensembl)

  df <- mapGene(refseq,to="ensembl_transcript")
  tx <- df[,2]
  expect_match(tx,"^ENST[0-9]+$")
  expect_false(any(is.na(tx)))
  expect_false("" %in% tx)

  df <- mapGene(refseq,to="entrez")
  entrez <- as.character(df[,2])
  expect_match(entrez,"^[0-9]+$")
  expect_false(any(is.na(entrez)))
  expect_false("" %in% entrez)

  df <- mapGene(refseq,to="hgnc_symbol")
  hgnc <- df[,2]
  expect_true(all(hgnc %in% gene_id_converter_hs$hgnc_symbol))
  expect_false(any(is.na(hgnc)))
  expect_false("" %in% hgnc)
})


test_that("mapGene() doesn't convert if the target gene id is same as current gene id",{
  ensembl <- rownames(exprs_dmap)
  df <- mapGene(ensembl,to="ensembl_gene")
  expect_equal(df[,1],df[,2])
})


test_that("mapGene() returns NULL if gene ids are not recognized",{
  expect_null(mapGene(c("123","ABC","EMM123")))
})
