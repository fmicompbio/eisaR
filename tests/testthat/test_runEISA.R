context("runEISA runs")
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
    res3 <- runEISA(cntEx, cntIn, cond, recalcLibSizeAfterFilt = TRUE, sizeFactor = "intron")
    expect_is(res0, "list")
    expect_is(res0se, "list")
    expect_equal(res0$contrasts, res0se$contrasts)
    expect_is(res1, "list")
    expect_is(res1se, "list")
    expect_equal(res1, res1se)
    expect_is(res2, "list")
    expect_is(res3, "list")
    expect_length(res1, 8L)
    expect_true(all(rownames(res1$DGEList) %in% rownames(cntEx)))
    ids <- intersect(rownames(res1$DGEList), rownames(res2$DGEList))
    expect_gt(cor(res1$contrasts[ids,"Dex"], res2$contrasts[ids,"Dex"]), 0.99)
    expect_gt(cor(res1$contrasts[ids,"Din"], res2$contrasts[ids,"Din"]), 0.99)
    expect_gt(cor(res1$contrasts[ids,"Dex.Din"], res2$contrasts[ids,"Dex.Din"]), 0.97)
    expect_gt(cor(res2$contrasts[ids,"Dex"], res3$contrasts[ids,"Dex"]), 0.99)
    
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

context("runEISA gives expected results")
test_that("runEISA() gives expected results", {
    # construct artificial example
    fracIntron <- 0.2
    ngenes <- 3000
    nsig <- 80
    nsig2 <- 20
    set.seed(5)

    # ... 2 conditions (ES, TN),
    #     2 replicates (exon counts, constant library size, just sampling noise)
    cond <- factor(c("ES","ES","TN","TN"))
    tmp <- readRDS(system.file("extdata", "Fig3abc_GSE33252_rawcounts_exonic.rds",
                               package = "eisaR"))[,-1]
    cntEx <- cbind(rmultinom(2, round(sum(tmp) / 4), rowMeans(tmp[, 1:2])),
                   rmultinom(2, round(sum(tmp) / 4), rowMeans(tmp[, 3:4])))
    colnames(cntEx) <- colnames(tmp)

    # ... select ngenes expressed genes
    selgenes <- sample(x = which(rowMeans(log2(cntEx + 1)) > 5.0),
                       size = ngenes, replace = FALSE)
    cntEx <- cntEx[selgenes, ]

    # ... intron counts as sub-sampled exon counts with fracIntron total counts
    nEx <- colSums(cntEx)
    cntIn <- cbind(rmultinom(2, round(sum(nEx[1:2]) / 2 * fracIntron), rowMeans(cntEx[, 1:2])),
                   rmultinom(2, round(sum(nEx[3:4]) / 2 * fracIntron), rowMeans(cntEx[, 3:4])))
    dimnames(cntIn) <- dimnames(cntEx)

    # ... nsig genes with significant a interaction
    #     - multiplicative effect, in order to keep the n_exon/N_total
    #       and n_intron/N_total ratios the same
    #     - half applied to ES and NP, each (in order to get both directions)
    selsig <- rownames(cntEx)[sample(x = nrow(cntEx), size = nsig, replace = FALSE)]
    lfcsig <- sample(x = 2:5, size = nsig, replace = TRUE)
    names(lfcsig) <- selsig
    i <- seq(1, round(nsig / 2))
    cntEx[selsig[i],   1:2] <- round(cntEx[selsig[i],   1:2] * 2^lfcsig[i])
    cntEx[selsig[-i],  3:4] <- round(cntEx[selsig[-i],  3:4] * 2^lfcsig[-i])

    # ... nsig2 genes with a sample specific effect
    #     - re-adjust counts of one replicate only, both for exons and introns
    #       (keeps the effect constant, but gives a sample-specific offset)
    selsig2 <- sample(x = selsig, size = nsig2, replace = FALSE)
    cntEx[selsig2,  4] <- round(cntEx[selsig2,  4] / 2^(lfcsig[match(selsig2, selsig)] / 2))
    cntIn[selsig2,  4] <- round(cntIn[selsig2,  4] / 2^(lfcsig[match(selsig2, selsig)] / 2))

    # run EISA (use gene filtering that is independent of model)
    res1 <- runEISA(cntEx, cntIn, cond, geneSelection = "Gaidatzis2015",
                    pscnt = 8, sizeFactor = "individual", modelSamples = FALSE)
    res2 <- runEISA(cntEx, cntIn, cond, geneSelection = "Gaidatzis2015",
                    pscnt = 8, sizeFactor = "individual", modelSamples = TRUE)

    # account for filtered genes
    ngenes <- nrow(res1$tab.ExIn)
    selsig <- intersect(selsig, rownames(res1$tab.ExIn)); nsig <- length(selsig)
    selsig2 <- intersect(selsig2, selsig); nsig2 <- length(selsig2)
    lfcsig <- lfcsig[selsig]

    # look at results
    # plot(res1$tab.ExIn[, "FDR"], res2$tab.ExIn[, "FDR"], log = "xy",
    #      col = (rownames(res1$tab.ExIn) %in% selsig) + (rownames(res1$tab.ExIn) %in% selsig2) + 1)
    # abline(a=0, b=1)
    # plot(res1$contrasts[selsig, "Dex.Din"], res2$contrasts[selsig, "Dex.Din"],
    #      col = (selsig %in% selsig2) + 1)
    # abline(a=0, b=1)
    # plotEISA(res1)
    # plotEISA(res2)
    t1 <- table(sample = rownames(res1$tab.ExIn) %in% selsig2,
                true.sig = rownames(res1$tab.ExIn) %in% selsig,
                test.sig = res1$tab.ExIn$FDR < 0.01) # close to 1/nsig -> one false positive
    t2 <- table(sample = rownames(res2$tab.ExIn) %in% selsig2,
                true.sig = rownames(res2$tab.ExIn) %in% selsig,
                test.sig = res2$tab.ExIn$FDR < 0.01) # close to 1/nsig -> one false positive

    # no false positives or false negatives (except for selsig2 genes in t1)
    nmissed1 <- t1["TRUE", "TRUE", "FALSE"] # expect to miss some genes with sample-specific effects
    expect_identical(as.vector(t1),
                     c(ngenes - nsig, 0L, 0L, nmissed1, 0L, 0L, nsig - nsig2, nsig2 - nmissed1))
    expect_identical(as.vector(t2),
                     c(ngenes - nsig, 0L, 0L, 0L, 0L, 0L, nsig - nsig2, nsig2))
})
