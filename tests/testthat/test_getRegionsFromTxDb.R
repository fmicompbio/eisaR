test_that("getRegionsFromTxDb() runs", {
    txdb <- AnnotationDbi::loadDb(system.file("extdata", "hg19sub.sqlite", package = "eisa"))
    regL <- getRegionsFromTxDb(txdb)

    # arguments
    expect_error(getRegionsFromTxDb(1:10))
    expect_error(getRegionsFromTxDb(txdb, exonExt = "a"))
    expect_error(getRegionsFromTxDb(txdb, strandedData = "a"))

    # expected results
    expect_true(is.list(regL))
    expect_length(regL, 2L)
    expect_equivalent(lengths(regL), c(20L, 4L))
    expect_true(all(countOverlaps(regL$exons) == 1L))
})
