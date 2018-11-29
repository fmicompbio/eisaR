test_that("plotEISA() runs", {
    cntEx <- readRDS(system.file("extdata", "Fig3abc_GSE33252_rawcounts_exonic.rds", package = "eisaR"))[,-1]
    cntIn <- readRDS(system.file("extdata", "Fig3abc_GSE33252_rawcounts_intronic.rds", package = "eisaR"))[,-1]
    cond <- factor(c("ES","ES","TN","TN"))
    res1 <- runEISA(cntEx, cntIn, cond, method = "published")

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_null(plotEISA(res1))
    expect_null(plotEISA(res1, contrast = "none"))

    dev.off()
    unlink(tf)
})
