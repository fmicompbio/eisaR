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
#'   that groups the samples (columns of \code{cntEx} and \code{cntIn} into two
#'   conditions. The contrast will be defined as secondLevel - firstLevel.
#' @param method One of \code{"published"} or \code{"new"}. If
#'   \code{"published"} (the default), sample normalization, gene filtering and
#'   calculation of contrasts is performed as described in Gaidatzis et al.
#'   2015, and the statistical analysis is based on \code{\link[edgeR]{glmFit}}
#'   and \code{\link[edgeR]{glmLRT}}.
#'   If \code{method="new"}, normalization will be performed using TMM as
#'   implemented in \code{\link[edgeR]{calcNormFactors}}, the contrast are
#'   calculated using \code{\link[edgeR]{predFC}}, and the statistical analysis
#'   will use the quasi-likelihood framework implemented in
#'   \code{\link[edgeR]{glmQLFit}} and \code{\link[edgeR]{glmQLFTest}}. The
#'   latter is often less stringent when selecting quantifyable genes, but more
#'   stringent wenn calling significant changes (especially with low numbers of
#'   replicates).
#' @param pscnt \code{numeric(1)} with pseudocount to add to read counts
#'   (default: 8). It is added to scaled read counts for
#'   \code{method="published"}, or used in \code{cpm(..., prior.count = pscnt)}
#'   and \code{predFC(..., prior.count = pscnt)} for \code{method="new"}.
#' @param ... additional arguments passed to the \code{\link[edgeR]{DGEList}}
#'   constructor,
#'   such as \code{lib.size} or \code{genes}.
#'
#' @return a \code{list} with elements \describe{
#'   \item{fracIn}{fraction intronic counts in each sample}
#'   \item{contrastName}{contrast name}
#'   \item{contrasts}{contrast matrix for quantifyable genes, with average log2
#'     fold-changes in exons (\code{Dex}), in introns (\code{Din}), and average
#'     difference between log2 fold-changes in exons and introns (\code{Dex.Din})}
#'   \item{DGEList}{\code{\link[edgeR]{DGEList}} object used in model fitting}
#'   \item{tab.cond}{statistical results for differential expression between
#'   conditions, based on both exonic and intronic counts}
#'   \item{tab.ExIn}{statisical results for differential changes between exonic
#'   and intronic contrast, an indication for post-transcriptional regulation.}
#'   \item{method}{the method that was used to run EISA}
#'   \item{pscnt}{the pseudocount that was used to run EISA}
#' }
#'
#' @references Analysis of intronic and exonic reads in RNA-seq data characterizes
#'   transcriptional and post-transcriptional regulation.
#'   Dimos Gaidatzis, Lukas Burger, Maria Florescu and Michael B. Stadler
#'   Nature Biotechnology, 2015 Jul;33(7):722-9. doi: 10.1038/nbt.3269.
#'
#' @seealso \code{\link[edgeR]{DGEList}} for \code{DGEList} object construction,
#'   \code{\link[edgeR]{calcNormFactors}} for normalization,
#'   \code{\link[edgeR]{glmFit}} for statisical analysis under \code{method = "published"}
#'   and \code{\link[edgeR]{glmQLFit}} statistical analysis under \code{method = "new"}.
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
#' @importFrom stats model.matrix
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay SummarizedExperiment
#'
#' @export
runEISA <- function(cntEx, cntIn, cond, method = c("published", "new"), 
                    pscnt = 8, ...) {
    # check arguments
    # ... count matrices
    if (is(cntEx, "SummarizedExperiment")) {
        cntIn <- SummarizedExperiment::assay(cntEx, "intron")
        cntEx <- SummarizedExperiment::assay(cntEx, "exon")
    }
    if (is.data.frame(cntEx))
        cntEx <- as.matrix(cntEx)
    stopifnot(is.matrix(cntEx))
    if (is.data.frame(cntIn))
        cntIn <- as.matrix(cntIn)
    stopifnot(is.matrix(cntIn))
    # ... consistency between cntEx and cntIn
    stopifnot(all(dim(cntEx) == dim(cntIn)))
    if (is.null(rownames(cntEx)))
        rownames(cntEx) <- as.character(seq.int(nrow(cntEx)))
    if (is.null(colnames(cntEx)))
        colnames(cntEx) <- as.character(seq.int(ncol(cntEx)))
    if (is.null(rownames(cntIn)))
        rownames(cntIn) <- as.character(seq.int(nrow(cntIn)))
    if (is.null(colnames(cntIn)))
        colnames(cntIn) <- as.character(seq.int(ncol(cntIn)))
    stopifnot(identical(dimnames(cntEx), dimnames(cntIn)))
    # ... conditions
    if (is.numeric(cond) || is.character(cond))
        cond <- factor(cond, levels = unique(cond))
    stopifnot(is.factor(cond))
    stopifnot(nlevels(cond) == 2L)
    stopifnot(length(cond) == ncol(cntEx))
    # ... method
    method <- match.arg(method)
    # ... pscnt
    stopifnot(is.numeric(pscnt) && length(pscnt) == 1)

    # fraction intronic
    fracIn <- colSums(cntIn) / (colSums(cntEx) + colSums(cntIn))

    # create DGEList
    cnt <- data.frame(Ex = cntEx, In = cntIn)
    y <- edgeR::DGEList(counts = cnt, ...)
    y <- edgeR::calcNormFactors(y)

    # identify quantifyable genes
    message("filtering quantifyable genes...", appendLF = FALSE)
    if (method == "published") {
        # scale counts to the mean library size separately for exons and introns
        Nex <- t(t(cntEx) / colSums(cntEx) * mean(colSums(cntEx)))
        Nin <- t(t(cntIn) / colSums(cntIn) * mean(colSums(cntIn)))

        # log transform (add pseudocount)
        if (pscnt != 8)
            warning("Using a 'pscnt' different from 8 deviates from method='published'")
        NLex <- log2(Nex + pscnt)
        NLin <- log2(Nin + pscnt)

        # Identify quantifyable genes
        quantGenes <- rownames(cntEx)[ rowMeans(NLex) > 5.0 & rowMeans(NLin) > 5.0 ]

    } else if (method == "new") {
        # Identify quantifyable genes (at least 5.0 reads per million in at least two samples)
        quantGenes <- rownames(cntEx)[ rowSums(edgeR::cpm(y, prior.count = pscnt) > 5.0) >= 2 ]
    }
    message("keeping ", length(quantGenes), " from ", nrow(y), " (",
            round(length(quantGenes) * 100 / nrow(y), 1), "%)")
    y <- y[quantGenes, ]
    y <- edgeR::calcNormFactors(y)

    # statistical analysis
    if (any(table(cond) < 2)) {
        warning("Need at least two replicates per condition to perform ",
                "statistical analysis. tt.cond and tt.ExIn will be empty.")
        tt.cond <- list(table = data.frame())
        tt.ExIn <- list(table = data.frame())
    } else {
        message("fitting statistical model...", appendLF = FALSE)
        cond2 <- rep(cond, 2L)
        region <- factor(rep(c("ex", "in"), each = ncol(cntEx)),
                         levels = c("in", "ex"))
        design <- model.matrix(~ region * cond2) # design matrix with interaction term
        rownames(design) <- colnames(cnt)
        y <- edgeR::estimateDisp(y, design)
        if (method == "published") {
            # use glmFit / glmLRT for method = "published"
            fit <- edgeR::glmFit(y, design)
            tst.cond <- edgeR::glmLRT(fit, contrast = c(0, 0, 1, 0.5))
            tst.ExIn <- edgeR::glmLRT(fit, coef = 4L)
        } else if (method == "new") {
            # use glmQLFit / glmQLFTest for method = "new"
            fit <- edgeR::glmQLFit(y, design)
            tst.cond <- edgeR::glmQLFTest(fit, contrast = c(0, 0, 1, 0.5))
            tst.ExIn <- edgeR::glmQLFTest(fit, coef = 4L)
        }
        tt.cond <- edgeR::topTags(tst.cond, n = nrow(y), sort.by = "none")
        tt.ExIn <- edgeR::topTags(tst.ExIn, n = nrow(y), sort.by = "none")
        message("done")
    }

    # calculate contrasts
    message("calculating contrasts...", appendLF = FALSE)
    contrastName <- paste(levels(cond)[2], "-", levels(cond)[1])
    if (method == "published") {
        i1 <- which(cond == levels(cond)[1])
        i2 <- which(cond == levels(cond)[2])
        Dex <- rowMeans(NLex[quantGenes, i2, drop = FALSE]) - 
            rowMeans(NLex[quantGenes, i1, drop = FALSE])
        Din <- rowMeans(NLin[quantGenes, i2, drop = FALSE]) - 
            rowMeans(NLin[quantGenes, i1, drop = FALSE])
        Dex.Din <- Dex - Din
    } else if (method == "new") {
        if (is.null(y$common.dispersion))
            stop("method='new' requires that dispersions can be calculated")
        lfc <- edgeR::predFC(y, design, prior.count = pscnt)
        rownames(lfc) <- rownames(y)
        Dex <- rowSums(lfc[, c(3, 4)])
        Din <- lfc[, 3]
        Dex.Din <- lfc[, 4]
    }
    message("done")

    ## return results
    return(list(fracIn = fracIn,
                contrastName = contrastName,
                contrasts = cbind(Dex = Dex, Din = Din, Dex.Din = Dex.Din),
                DGEList = y, tab.cond = tt.cond$table, tab.ExIn = tt.ExIn$table,
                method = method, pscnt = pscnt))
}

