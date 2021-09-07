test_that("getRegionsFromTxDb() runs", {
    requireNamespace("GenomicFeatures")
    txdb <- AnnotationDbi::loadDb(system.file("extdata", "hg19sub.sqlite", package = "eisaR"))
    regL <- getRegionsFromTxDb(txdb)

    # arguments
    expect_error(getRegionsFromTxDb(1:10))
    expect_error(getRegionsFromTxDb(txdb, exonExt = "a"))
    expect_error(getRegionsFromTxDb(txdb, strandedData = "a"))
    
    # required installed packages
    # ... set new lib paths
    old <- .libPaths()
    td <- tempfile(pattern = "Rlib")
    dir.create(td)
    .libPaths(c(td, old[length(old)]), include.site = FALSE)
    unloadNamespace("ensembldb")
    unloadNamespace("GenomicFeatures")
    # ... test
    expect_error(getRegionsFromTxDb(structure("dummy", class = "TxDb")))
    # ... clean up
    unlink(td, recursive = TRUE, force = TRUE)
    .libPaths(old)
    requireNamespace("GenomicFeatures")
    requireNamespace("ensembldb")
    
    # expected results
    expect_true(is.list(regL))
    expect_length(regL, 2L)
    expect_equivalent(lengths(regL), c(20L, 4L))
    expect_true(all(countOverlaps(regL$exons) == 1L))
    
    # same for EnsDb object
    # ... load GRanges with annotation (Y)
    load(system.file("YGRanges.RData", package = "ensembldb"))
    # ... create EnsDb
    suppressWarnings(DB <- ensembldb::ensDbFromGRanges(Y, path = tempdir(), version = 75,
                                                       organism = "Homo_sapiens"))
    edb <- ensembldb::EnsDb(DB)
    # ... create corresponding TxDb
    suppressWarnings(txdb <- GenomicFeatures::makeTxDbFromGRanges(Y))
    # ... compare annotation
    expect_identical(length(GenomicFeatures::genes(edb)),
                     length(GenomicFeatures::genes(txdb)))
    expect_identical(length(GenomicFeatures::exons(edb)),
                     length(GenomicFeatures::exons(txdb)))
    # ... compare results
    regL1 <- getRegionsFromTxDb(edb)
    regL2 <- getRegionsFromTxDb(txdb)
    expect_identical(regL1, regL2)
})
