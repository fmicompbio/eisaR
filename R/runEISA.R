#' @title Run Exon-Intron Split Analysis.
#'
#' @description Starting from count tables with exonic and intronic counts
#'   for two conditions, perform all the steps in EISA (normalize, identify
#'   quantifyable genes, calculate contrasts and their significance).
#'
#' @author Michael Stadler
#'
#' @param cntEx Gene by sample \code{matrix} with exonic counts, OR a
#'   \code{SummarizedExperiment} with two assays named \code{exon} and
#'   \code{intron}, containing exonic and intronic counts, respectively. If
#'   \code{cntEx} is a \code{SummarizedExperiment}, \code{cntIn} will be
#'   disregarded.
#' @param cntIn Gene by sample \code{matrix} with intronic counts. Must have the
#'   same structure as \code{cntEx} (same number and order of rows and columns)
#'   if \code{cntEx} is a matrix. Will be disregarded if \code{cntEx} is a
#'   \code{SummarizedExperiment}.
#' @param cond \code{numeric}, \code{character} or \code{factor} with two levels
#'   that groups the samples (columns of \code{cntEx} and \code{cntIn}) into two
#'   conditions. The contrast will be defined as secondLevel - firstLevel.
#' @param method One of \code{NULL} (the default) or \code{"Gaidatzis2015"}. If
#'   \code{"Gaidatzis2015"}, gene filtering, statistical analysis and
#'   calculation of contrasts is performed as described in Gaidatzis et al.
#'   2015, and the statistical analysis is based on \code{\link[edgeR]{glmFit}}
#'   and \code{\link[edgeR]{glmLRT}}. This is done by setting the arguments
#'   \code{modelSamples}, \code{geneSelection}, \code{effects} and \code{statFramework}
#'   to appropriate values (see details), overriding the defaults or any value
#'   passed to these arguments. If \code{NULL}, the default values of the arguments
#'   will be used instead (recommended).
#' @param modelSamples Whether to include a sample identifier in the design matrix
#'   of the statistical model. If \code{TRUE}, potential sample effects
#'   that affect both exonic and intronic counts of that sample will be taken
#'   into account, which could result in higher sensitivity (default: FALSE).
#' @param geneSelection Controls how to select quantifyable genes. One of the
#'   following:\describe{
#'       \item{\code{"filterByExpr"}: }{(default) First, counts are normalized using
#'       \code{\link[edgeR]{calcNormFactors}}, treating intonic and exonic counts
#'       as individual samples. Then, \code{\link[edgeR]{filterByExpr}} is used
#'       with default parameters to select quantifyable genes.}
#'       \item{\code{"none"}: }{This will use all the genes provided in the count
#'       tables, assuming that an appropriate selection of quantifyable genes has
#'       already been done.}
#'       \item{\code{"Gaidatzis2015"}: }{First, intronic and exonic counts are
#'       linearly scaled to the mean library size (estimated as the sum of all
#'       intronic or exonic counts, respectively). Then, quantifyable genes are
#'       selected as the genes with counts \code{x} that fulfill
#'       \code{log2(x + 8) > 5} in both exons and introns.}
#'   }
#' @param statFramework Selects the framework within \code{edgeR} that is used
#'   for the statistical analysis. One of:\describe{
#'       \item{\code{"QLF"}: }{(default) Quasi-likelihood F-test using
#'       \code{\link[edgeR]{glmQLFit}} and \code{\link[edgeR]{glmQLFTest}}. This
#'       framework is highly recommended as it gives stricter error rate control
#'       by accounting for the uncertainty in dispersion estimation.}
#'       \item{\code{"LRT"}: }{Likelyhood ratio test using \code{\link[edgeR]{glmFit}}
#'       and \code{\link[edgeR]{glmLRT}}}.
#'   }
#' @param effects How the effects (contrasts or log2 fold-changes) are calculated.
#'   One of:\describe{
#'       \item{\code{"predFC"}: }{(default) Fold-changes are calculated using
#'       the fitted model with \code{\link[edgeR]{predFC}} with
#'       \code{prior.count = pscnt}. Please note that if a sample factor is
#'       included in the model (\code{modelSamples=TRUE}), effects cannot be
#'       obtained from that model. In that case, effects are obtained from a
#'       simpler model without sample effects.}
#'       \item{\code{"Gaidatzis2015"}: }{Fold-changes are calculated using the
#'       formula \code{log2((x + pscnt)/(y + pscnt))}. If \code{pscnt} is not
#'       set to 8, \code{runEISA} will warn that this deviates from the method
#'       used in Gaidatzis et al., 2015.}
#'   }
#' @param pscnt \code{numeric(1)} with pseudocount to add to read counts
#'   (default: 2). For \code{method = "Gaidatzis2015"}, it is set to 8.
#'   It is added to scaled read counts used in \code{geneSelection = "Gaidatzis2015"}
#'   and \code{effects = "Gaidatzis2015"}, or else used in \code{cpm(..., prior.count = pscnt)}
#'   and \code{predFC(..., prior.count = pscnt)}.
#' @param ... additional arguments passed to the \code{\link[edgeR]{DGEList}}
#'   constructor, such as \code{lib.size} or \code{genes}.
#'
#' @details Setting \code{method = "Gaidatzis2015"} has precedence over other
#'   argument values and corresponds to setting:
#'   \code{modelSamples = FALSE, geneSelection = "Gaidatzis2015",
#'   statFramework = "LRT", effects = "Gaidatzis2015", pscnt = 8}.
#'   
#' @return a \code{list} with elements \describe{
#'   \item{fracIn}{fraction intronic counts in each sample}
#'   \item{contrastName}{contrast name}
#'   \item{contrasts}{contrast matrix for quantifyable genes, with average log2
#'     fold-changes in exons (\code{Dex}), in introns (\code{Din}), and average
#'     difference between log2 fold-changes in exons and introns (\code{Dex.Din})}
#'   \item{DGEList}{\code{\link[edgeR]{DGEList}} object used in model fitting}
#'   \item{tab.ExIn}{statisical results for differential changes between exonic
#'   and intronic contrast, an indication for post-transcriptional regulation.}
#'   \item{params}{a \code{list} with parameter values used to run EISA}
#' }
#'
#' @references Analysis of intronic and exonic reads in RNA-seq data characterizes
#'   transcriptional and post-transcriptional regulation.
#'   Dimos Gaidatzis, Lukas Burger, Maria Florescu and Michael B. Stadler
#'   Nature Biotechnology, 2015 Jul;33(7):722-9. doi: 10.1038/nbt.3269.
#'
#' @seealso \code{\link[edgeR]{DGEList}} for \code{DGEList} object construction,
#'   \code{\link[edgeR]{calcNormFactors}} for normalization,
#'   \code{\link[edgeR]{filterByExpr}} for gene selection,
#'   \code{\link[edgeR]{glmFit}} and \code{\link[edgeR]{glmQLFit}} for statistical
#'   analysis.
#'
#' @examples
#' cntEx <- readRDS(system.file("extdata", "Fig3abc_GSE33252_rawcounts_exonic.rds",
#'                              package = "eisaR"))[,-1]
#' cntIn <- readRDS(system.file("extdata", "Fig3abc_GSE33252_rawcounts_intronic.rds",
#'                              package = "eisaR"))[,-1]
#' cond <- factor(c("ES","ES","TN","TN"))
#' res <- runEISA(cntEx, cntIn, cond)
#' plotEISA(res)
#'
#' @import edgeR
#' @importFrom limma nonEstimable
#' @importFrom stats model.matrix
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay assayNames SummarizedExperiment
#'
#' @export
runEISA <- function(cntEx, cntIn, cond, method = NULL, 
                    modelSamples = FALSE,
                    geneSelection = c("filterByExpr", "none", "Gaidatzis2015"),
                    statFramework = c("QLF", "LRT"),
                    effects = c("predFC", "Gaidatzis2015"),
                    pscnt = 2, ...) {
    # check arguments
    # ... count matrices
    if (is(cntEx, "SummarizedExperiment")) {
        if (all(c("exon","intron") %in% SummarizedExperiment::assayNames(cntEx))) {
            cntIn <- SummarizedExperiment::assay(cntEx, "intron")
            cntEx <- SummarizedExperiment::assay(cntEx, "exon")
        } else if (all(c("spliced","unspliced") %in% SummarizedExperiment::assayNames(cntEx))) {
            cntIn <- SummarizedExperiment::assay(cntEx, "unspliced")
            cntEx <- SummarizedExperiment::assay(cntEx, "spliced")
        } else {
            stop("'cntEx' needs to have assayNames 'intron'/'exon' or 'unspliced'/'spliced'.")
        }
    }
    if (is.data.frame(cntEx))
        cntEx <- as.matrix(cntEx)
    if (is.data.frame(cntIn))
        cntIn <- as.matrix(cntIn)
    stopifnot(exprs = {
        is.matrix(cntEx)
        is.matrix(cntIn)
    })
    # ... consistency between cntEx and cntIn
    stopifnot(all(dim(cntEx) == dim(cntIn)))
    nsmpls <- ncol(cntEx)
    if (is.null(rownames(cntEx)))
        rownames(cntEx) <- as.character(seq.int(nrow(cntEx)))
    if (is.null(colnames(cntEx)))
        colnames(cntEx) <- as.character(seq.int(nsmpls))
    if (is.null(rownames(cntIn)))
        rownames(cntIn) <- as.character(seq.int(nrow(cntIn)))
    if (is.null(colnames(cntIn)))
        colnames(cntIn) <- as.character(seq.int(ncol(cntIn)))
    stopifnot(identical(dimnames(cntEx), dimnames(cntIn)))
    # ... conditions
    if (is.numeric(cond) || is.character(cond))
        cond <- factor(cond, levels = unique(cond))
    # ... valid arguments
    geneSelection <- match.arg(geneSelection)
    statFramework <- match.arg(statFramework)
    effects <- match.arg(effects)
    stopifnot(exprs = {
        # cond
        is.factor(cond)
        nlevels(cond) == 2L
        length(cond) == nsmpls
        # method
        is.null(method) || method %in% c("Gaidatzis2015")
        # modelSamples
        is.logical(modelSamples)
        length(modelSamples) == 1L
        # pscnt
        is.numeric(pscnt)
        length(pscnt) == 1L
    })

    # override arguments for Gaidatzis2015
    if (!is.null(method) && method == "Gaidatzis2015") {
        message("setting parameters according to Gaidatzis et al., 2015")
        modelSamples <- FALSE
        geneSelection <- "Gaidatzis2015"
        statFramework <- "LRT"
        effects <- "Gaidatzis2015"
        pscnt <- 8
    }
    
    # fraction intronic
    fracIn <- colSums(cntIn) / (colSums(cntEx) + colSums(cntIn))

    # create DGEList
    cnt <- data.frame(Ex = cntEx, In = cntIn)
    y <- edgeR::DGEList(counts = cnt, ...)
    y <- edgeR::calcNormFactors(y)

    # create design matrix
    cond2 <- rep(cond, 2L)
    region <- factor(rep(c("ex", "in"), each = nsmpls),
                     levels = c("in", "ex"))
    smpl <- factor(rep(sprintf("s%03d", seq.int(nsmpls)), 2))
    if (modelSamples) {
        dsgn <- model.matrix(~ smpl + region * cond2)
        # need to remove a coefficient to make the design full rank
        toRemove <- limma::nonEstimable(dsgn)
        dsgn <- dsgn[, -match(toRemove, colnames(dsgn))]
    } else {
        dsgn <- model.matrix(~ region * cond2)
    }
    rownames(dsgn) <- colnames(cnt)
    
    # identify quantifyable genes
    if (geneSelection == "none") {
        message("skip filtering for quantifyable genes...", appendLF = FALSE)
        NLex <- edgeR::cpm(y[, seq.int(nsmpls)], log = TRUE, prior.count = pscnt)
        NLin <- edgeR::cpm(y[, nsmpls + seq.int(nsmpls)], log = TRUE, prior.count = pscnt)
        
    } else {
        message("filtering quantifyable genes...", appendLF = FALSE)

        if (geneSelection == "filterByExpr") {
            quantGenes <- rownames(cntEx)[ edgeR::filterByExpr(y[, seq.int(nsmpls)],
                                                               design = dsgn[seq.int(nsmpls), ]) &
                                           edgeR::filterByExpr(y[, nsmpls + seq.int(nsmpls)],
                                                               design = dsgn[nsmpls + seq.int(nsmpls), ]) ]
            NLex <- edgeR::cpm(y[, seq.int(nsmpls)], log = TRUE, prior.count = pscnt)
            NLin <- edgeR::cpm(y[, nsmpls + seq.int(nsmpls)], log = TRUE, prior.count = pscnt)
            
        } else if (geneSelection == "Gaidatzis2015") {
            # scale counts to the mean library size separately for exons and introns
            Nex <- t(t(cntEx) / colSums(cntEx) * mean(colSums(cntEx)))
            Nin <- t(t(cntIn) / colSums(cntIn) * mean(colSums(cntIn)))

            # log transform (add pseudocount)
            if (pscnt != 8)
                warning("Using a 'pscnt' different from 8 deviates from geneSelection='Gaidatzis2015'")
            NLex <- log2(Nex + pscnt)
            NLin <- log2(Nin + pscnt)

            # Identify quantifyable genes
            quantGenes <- rownames(cntEx)[ rowMeans(NLex) > 5.0 & rowMeans(NLin) > 5.0 ]

        }
        message("keeping ", length(quantGenes), " from ", nrow(y), " (",
                round(length(quantGenes) * 100 / nrow(y), 1), "%)")
        y <- y[quantGenes, ]
        y <- edgeR::calcNormFactors(y)
        NLex <- NLex[quantGenes, ]
        NLin <- NLin[quantGenes, ]
    }

    # statistical analysis
    if (any(table(cond) < 2)) {
        warning("Need at least two replicates per condition to perform ",
                "statistical analysis. 'ExIn' result will be empty.")
        tt.ExIn <- list(table = data.frame())

    } else {
        message("fitting statistical model...", appendLF = FALSE)
        y <- edgeR::estimateDisp(y, dsgn)
        if (statFramework == "QLF") {
            fit <- edgeR::glmQLFit(y, dsgn)
            tst.ExIn <- edgeR::glmQLFTest(fit, coef = ncol(dsgn))
        } else if (statFramework == "LRT") {
            fit <- edgeR::glmFit(y, dsgn)
            tst.ExIn <- edgeR::glmLRT(fit, coef = ncol(dsgn))
        }
        tt.ExIn <- edgeR::topTags(tst.ExIn, n = nrow(y), sort.by = "none")
        message("done")
    }

    # calculate contrasts
    message("calculating contrasts...", appendLF = FALSE)
    contrastName <- paste(levels(cond)[2], "-", levels(cond)[1])
    if (effects == "predFC") {
        if (is.null(y$common.dispersion))
            stop("effects='predFC' requires a fitted model - rerun with effects='Gaidatzis2015'")
        lfc <- edgeR::predFC(y, dsgn, prior.count = pscnt)
        if (modelSamples) {
            message("fitting model without sample factor...", appendLF = FALSE)
            dsgn2 <- model.matrix(~ region * cond2)
            y2 <- edgeR::estimateDisp(y, dsgn2)
            if (statFramework == "QLF") {
                fit2 <- edgeR::glmQLFit(y2, dsgn2)
            } else if (statFramework == "LRT") {
                fit2 <- edgeR::glmFit(y2, dsgn2)
            }
            lfc2 <- edgeR::predFC(y2, dsgn2, prior.count = pscnt)
        } else {
            lfc2 <- lfc
        }
        rownames(lfc) <- rownames(lfc2) <- rownames(y)
        Dex <- rowSums(lfc2[, c(3, 4)])
        Din <- lfc2[, 3]
        Dex.Din <- lfc[, ncol(lfc)]
        # remark: for modelSamples=TRUE, should the interaction effect be estimated...
        #         - from the simpler model (as the condition effects, for consistency among effects)
        #         - from the full model (for consistency with the interaction FDR and topTags table) -> current implementation
    } else if (effects == "Gaidatzis2015") {
        i1 <- which(cond == levels(cond)[1])
        i2 <- which(cond == levels(cond)[2])
        Dex <- rowMeans(NLex[quantGenes, i2, drop = FALSE]) - 
            rowMeans(NLex[quantGenes, i1, drop = FALSE])
        Din <- rowMeans(NLin[quantGenes, i2, drop = FALSE]) - 
            rowMeans(NLin[quantGenes, i1, drop = FALSE])
        Dex.Din <- Dex - Din
    }
    message("done")

    ## return results
    return(list(fracIn = fracIn,
                contrastName = contrastName,
                contrasts = cbind(Dex = Dex, Din = Din, Dex.Din = Dex.Din),
                DGEList = y,
                tab.ExIn = tt.ExIn$table,
                params = list(method = method, modelSamples = modelSamples,
                              geneSelection = geneSelection, statFramework = statFramework,
                              effects = effects, pscnt = pscnt)))
}

