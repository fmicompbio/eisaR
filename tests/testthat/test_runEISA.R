test_that("runEISA() runs", {
    cntEx <- readRDS(system.file("extdata", "Fig3abc_GSE33252_rawcounts_exonic.rds", package = "eisaR"))[,-1]
    cntIn <- readRDS(system.file("extdata", "Fig3abc_GSE33252_rawcounts_intronic.rds", package = "eisaR"))[,-1]
    cond <- factor(c("ES","ES","TN","TN"))
    cntSE1 <- SummarizedExperiment::SummarizedExperiment(assays = list(exon = cntEx, intron = cntIn))
    cntSE2 <- SummarizedExperiment::SummarizedExperiment(assays = list(spliced = cntEx, unspliced = cntIn))
    
    # arguments
    expect_error(runEISA(data.frame("a")))
    expect_error(runEISA(rbind(1:2, 3:4), "b"))
    expect_error(runEISA(rbind(1:2, 3:4), rbind(1:2, 3:4), c("a","b","c")))
    expect_warning(runEISA(rbind(1:2, 3:4), rbind(1:2, 3:4), c("a","b"), effects = "Gaidatzis2015"))
    expect_warning(runEISA(cntEx, cntIn, cond, geneSelection = "Gaidatzis2015", pscnt = 4))
    expect_error(runEISA(SummarizedExperiment(assays = list(exon = cntEx))))

    # expected results
    res0 <- runEISA(cntEx[1:1000,], cntIn[1:1000,], cond,
                    modelSamples = TRUE, geneSelection = "none", statFramework = "LRT")
    res0se <- runEISA(cntEx = cntSE2[1:1000,], cntIn = NULL, cond,
                      modelSamples = TRUE, geneSelection = "none")
    res1 <- runEISA(cntEx, cntIn, cond, method = "Gaidatzis2015")
    res1se <- runEISA(cntEx = cntSE1, cntIn = NULL, cond, method = "Gaidatzis2015")
    res2 <- runEISA(cntEx, cntIn, cond, method = NULL, modelSamples = FALSE, sizeFactor = "individual")
    expect_is(res0, "list")
    expect_is(res0se, "list")
    expect_equal(res0$contrasts, res0se$contrasts)
    expect_is(res1, "list")
    expect_is(res1se, "list")
    expect_equal(res1, res1se)
    expect_is(res2, "list")
    expect_length(res1, 8L)
    expect_true(all(rownames(res1$DGEList) %in% rownames(cntEx)))
    ids <- intersect(rownames(res1$DGEList), rownames(res2$DGEList))
    expect_gt(cor(res1$contrasts[ids,"Dex"], res2$contrasts[ids,"Dex"]), 0.99)
    expect_gt(cor(res1$contrasts[ids,"Din"], res2$contrasts[ids,"Din"]), 0.99)
    expect_gt(cor(res1$contrasts[ids,"Dex.Din"], res2$contrasts[ids,"Dex.Din"]), 0.97)

    # one replicate per condition
    expect_warning(res1 <- runEISA(cntEx[, c(1, 3)], cntIn[, c(1, 3)],
                                   cond[c(1, 3)], method = "Gaidatzis2015"))
    expect_is(res1, "list")
    expect_length(res1, 8L)
    expect_named(res1, c("fracIn", "contrastName", "contrasts", "DGEList",
                         "tab.ExIn", "contr.ExIn", "designMatrix", "params"))
    expect_is(res1$tab.ExIn, "data.frame")
    expect_equal(nrow(res1$tab.ExIn), 0)
    expect_error(plotEISA(res1))
    expect_error(suppressWarnings(runEISA(cntEx[, c(1, 3)], cntIn[, c(1, 3)],
                                          cond[c(1, 3)], method = NULL)))
})
