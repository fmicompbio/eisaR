#' @title Visualize the results from an exon-intron split analysis.
#'
#' @description \code{plotEISA} takes the return value from \code{\link{runEISA}}
#'   and generates a scatterplot of intronic versus exonic changes.
#'
#' @author Michael Stadler
#'
#' @param x \code{list} with EISA results, typically the return value from \code{\link{runEISA}}
#' @param contrast one of \code{"Dex.Din"}, \code{"cond"} or \code{"none"}. If \code{"Dex.Din"}
#'   (the default), genes that significantly differ between exonic and intronic changes
#'   are highlighted. If \code{contrast="cond"}, genes that differ between conditions
#'   are highlighted. \code{"none"} turns off gene highlighting.
#' @param minLfc \code{NULL} or \code{numeric(1)} with the minimal absolute log2
#'   fold change to color a gene. If \code{NULL} (the default), no fold changes
#'   are not used to select genes for highlighting.
#' @param maxFDR \code{numeric(1)} with maximal false discovery rate for gene
#'   highlighting.
#' @param genecolors Vector of length three specifying the colors to use for
#'   genes that are significantly up, down or unchanged.
#' @param ... further arguments past to \code{plot()}
#'
#' @return `NULL` (invisibly)
#'
#' @examples
#'   # see the help for runEISA() for a full example
#'
#' @importFrom graphics plot legend
#'
#' @export
plotEISA <- function(x, contrast = c("ExIn", "cond", "none"),
                     minLfc = NULL, maxFDR = 0.05,
                     genecolors = c("#E41A1CFF", "#4DAF4AFF", "#22222244"), ...) {
    # check arguments
    contrast <- match.arg(contrast)
    sigtab <- switch(contrast, ExIn = x$tab.ExIn, cond = x$tab.cond, none = data.frame())
    if (nrow(sigtab) == 0 && contrast != "none")
        stop("'x' does not contain the requested statistics and can only ",
             "be plotted using contrast = 'none'. Note that at least two ",
             "replicates per condition are required to run the statistical ",
             "testing.")
    if (is.null(minLfc))
        minLfc <- 0
    stopifnot(is.numeric(minLfc) && length(minLfc) == 1L)
    stopifnot(is.numeric(maxFDR) && length(maxFDR) == 1L && maxFDR >= 0 && maxFDR <= 1.0)
    stopifnot(is.character(genecolors) && length(genecolors) == 3L)

    # identify gene to highlight
    if (contrast == "none") {
        sig <- rep(FALSE, nrow(x$contrasts))
        sig.dir <- numeric(0)
    } else {
        sig <- abs(sigtab$logFC) >= minLfc & sigtab$FDR <= maxFDR
        sig.dir <- sign(sigtab$logFC[sig])
        message("identified ", sum(sig), " genes to highlight")
    }

    # set graphical parameters (user-defined colors take precedence)
    dotsL <- list(...)
    dotsL$x <- x$contrasts[, "Din"]
    dotsL$y <- x$contrasts[, "Dex"]
    if (!"pch" %in% names(dotsL))
        dotsL$pch <- 20
    if (!"col" %in% names(dotsL))
        dotsL$col <- ifelse(sig, ifelse(sigtab$logFC > 0, genecolors[1], genecolors[2]), genecolors[3])
    if (!"xlab" %in% names(dotsL))
        dotsL$xlab <- substitute(expression(paste(Delta, "intron (", cn, ")")), list(cn = x$contrastName))
    if (!"ylab" %in% names(dotsL))
        dotsL$ylab <- substitute(expression(paste(Delta, "exon (", cn, ")")), list(cn = x$contrastName))

    # Delta I vs. Delta E
    do.call(plot, dotsL)
    if (contrast != "none") {
        legend(x = "bottomright", bty = "n", pch = 20, col = genecolors[1:2],
               legend = sprintf("%s (%d)", c("Up","Down"), c(sum(sig.dir == 1), sum(sig.dir == -1))))
    }

    return(invisible(NULL))
}
