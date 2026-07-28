test_that("plotEISA() runs", {
    cntEx <- readRDS(system.file("extdata",
                                 "Fig3abc_GSE33252_rawcounts_exonic.rds",
                                 package = "eisaR"))[, -1L]
    cntIn <- readRDS(system.file("extdata",
                                 "Fig3abc_GSE33252_rawcounts_intronic.rds",
                                 package = "eisaR"))[, -1L]
    cond <- factor(c("ES", "ES", "TN", "TN"))
    res1 <- runEISA(cntEx, cntIn, cond, method = "Gaidatzis2015")
    expect_warning(res2 <- runEISA(cntEx[, c(1L, 3L)], cntIn[, c(1L, 3L)],
                                   cond[c(1L, 3L)], method = "Gaidatzis2015"))

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_true(ggplot2::is_ggplot(
        plotEISA(res1, genecolors = c("orange", "yellow", "purple"))))
    expect_true(ggplot2::is_ggplot(
        plotEISA(res1, contrast = "none")))
    expect_warning(plotEISA(res1, contrast = "none", cex = 2), 
                   "Ignoring unknown arguments cex")
    expect_error(plotEISA(res2))

    dev.off()
    unlink(tf)
})
