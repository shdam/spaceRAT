test_that("convertGeneName() converts between Ensembl_gene and Entrez",{
  exprs_entrez <- convertGeneName(exprs_dmap,to="entrez")
  expect_true(all(rownames(exprs_entrez) %in% gene_id_converter_hs$entrez))
  expect_true(ncol(exprs_entrez)==ncol(exprs_dmap))
  expect_true(nrow(exprs_entrez)<=nrow(exprs_dmap))
  expect_false( "" %in% rownames(exprs_entrez))

  exprs_ensembl <- convertGeneName(exprs_entrez,to="ensembl_gene")
  expect_match(rownames(exprs_ensembl),"^ENSG[0-9]+$")
  expect_true(ncol(exprs_ensembl)==ncol(exprs_entrez))
  expect_true(nrow(exprs_ensembl)<=nrow(exprs_entrez))
  expect_false( "" %in% rownames(exprs_ensembl))
})


test_that("convertGeneName() converts between Ensembl_gene and hgnc_symbol",{
  exprs_symbol <- convertGeneName(exprs_dmap,to="hgnc_symbol")
  expect_true(all(rownames(exprs_symbol) %in% gene_id_converter_hs$hgnc_symbol ))
  expect_true(ncol(exprs_symbol)==ncol(exprs_dmap))
  expect_true(nrow(exprs_symbol)<=nrow(exprs_dmap))
  expect_false( "" %in% rownames(exprs_symbol))

  exprs_ensembl <- convertGeneName(exprs_symbol,to="ensembl_gene")
  expect_match(rownames(exprs_ensembl),"^ENSG[0-9]+$")
  expect_true(ncol(exprs_ensembl)==ncol(exprs_symbol))
  expect_true(nrow(exprs_ensembl)<=nrow(exprs_symbol))
})


test_that("convertGeneName() converts from Ensembl_transcript to Ensembl_gene",{
  transcripts <- unique(gene_id_converter_hs$ensembl_transcript_id[!is.na(gene_id_converter_hs$ensembl_transcript_id)])[1:5000]
  dat <- matrix(1,nrow=5000,ncol=3)
  rownames(dat) <- transcripts
  df <- convertGeneName(dat)
  expect_true(ncol(df)==ncol(dat))
  expect_true(nrow(df)<=nrow(dat))
  expect_true(all(rownames(df) %in% gene_id_converter_hs$ensembl_gene))
})


test_that("convertGeneName() converts from refseq_mrna to Ensembl_gene",{
  transcripts <- unique(gene_id_converter_hs$refseq_mrna[!is.na(gene_id_converter_hs$refseq_mrna)])[1:5000]
  dat <- matrix(1,nrow=5000,ncol=3)
  rownames(dat) <- transcripts
  df <- convertGeneName(dat)
  expect_true(ncol(df)==ncol(dat))
  expect_true(nrow(df)<=nrow(dat))
  expect_true(all(rownames(df) %in% gene_id_converter_hs$ensembl_gene))
})
