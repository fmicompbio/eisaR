test_that("plotEISA() runs", {
    cntEx <- readRDS(system.file("extdata", "Fig3abc_GSE33252_rawcounts_exonic.rds", package = "eisaR"))[,-1]
    cntIn <- readRDS(system.file("extdata", "Fig3abc_GSE33252_rawcounts_intronic.rds", package = "eisaR"))[,-1]
    cond <- factor(c("ES","ES","TN","TN"))
    res1 <- runEISA(cntEx, cntIn, cond, method = "Gaidatzis2015")
    expect_warning(res2 <- runEISA(cntEx[, c(1, 3)], cntIn[, c(1, 3)], cond[c(1, 3)], method = "Gaidatzis2015"))

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_null(plotEISA(res1))
    expect_null(plotEISA(res1, contrast = "none"))
    expect_error(plotEISA(res2))

    dev.off()
    unlink(tf)
})
