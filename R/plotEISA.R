#' @title Visualize the results from an exon-intron split analysis.
#'
#' @description \code{plotEISA} takes the return value from \code{\link{runEISA}}
#'     and generates a scatterplot of intronic versus exonic changes.
#'
#' @author Michael Stadler, Charlotte Soneson
#'
#' @param x \code{list} with EISA results, typically the return value from 
#'     \code{\link{runEISA}}
#' @param contrast one of \code{"ExIn"} or \code{"none"}. If \code{"ExIn"}
#'     (the default), genes that significantly differ between exonic and 
#'     intronic changes are highlighted. \code{"none"} turns off gene 
#'     highlighting.
#' @param minLfc \code{NULL} or \code{numeric(1)} with the minimal absolute log2
#'     fold change to color a gene. If \code{NULL} (the default), no fold changes
#'     are not used to select genes for highlighting.
#' @param maxFDR \code{numeric(1)} with maximal false discovery rate for gene
#'     highlighting.
#' @param genecolors Vector of length three specifying the colors to use for
#'     genes that are significantly up, down or unchanged.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#'   # see the help for runEISA() for a full example
#'
#' @importFrom ggplot2 ggplot aes geom_abline geom_point coord_fixed 
#' @importFrom ggplot2 scale_color_manual labs theme_bw theme element_blank
#' @importFrom rlang .data
#'
#' @export
plotEISA <- function(x, contrast = c("ExIn", "none"),
                     minLfc = NULL, maxFDR = 0.05,
                     genecolors = c("#E41A1C", "#497AB3", "#222222")) {
    # check arguments
    contrast <- match.arg(contrast)
    contrastName <- ifelse("contrastName" %in% names(x),
                           paste0(" (", x$contrastName, ")"), "")
    sigtab <- switch(contrast, ExIn = x$tab.ExIn, none = data.frame())
    if (contrast != "none") {
        if (nrow(sigtab) == 0L) {
            stop("'x' does not contain the requested statistics and can only ",
                 "be plotted using contrast = 'none'. Note that at least two ",
                 "replicates per condition are required to run the statistical ",
                 "testing.")
        }
        if (nrow(sigtab) != nrow(x$contrasts)) {
            stop("x$tab.ExIn and x$contrasts have different numbers of rows.")
        }
    }
    if (is.null(minLfc)) {
        minLfc <- 0.0
    }
    stopifnot(exprs = {
        is.numeric(minLfc)
        length(minLfc) == 1L
        is.numeric(maxFDR)
        length(maxFDR) == 1L
        maxFDR >= 0.0
        maxFDR <= 1.0
        is.character(genecolors)
        length(genecolors) == 3L
    })

    # identify genes to highlight
    plotdf <- as.data.frame(x$contrasts)
    plotdf$sig <- ""
    if (contrast != "none") {
        plotdf$sig <- ifelse(
            sigtab$logFC >= minLfc & sigtab$FDR <= maxFDR, 
            paste0("Up (", sum(sigtab$logFC >= minLfc & sigtab$FDR <= maxFDR), ")"), 
            ifelse(sigtab$logFC <= -minLfc & sigtab$FDR <= maxFDR, 
                   paste0("Down (", sum(sigtab$logFC <= -minLfc & sigtab$FDR <= maxFDR), ")"),
                   ""))
        plotdf$sig <- factor(plotdf$sig, levels = rev(sort(unique(plotdf$sig))))
        message("identified ", sum(plotdf$sig != ""), " genes to highlight")
    }
    
    ggplot(plotdf, aes(x = .data$Din, y = .data$Dex)) + 
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                    color = "grey") + 
        geom_point(size = 0.1, color = genecolors[3L]) + 
        geom_point(data = plotdf[plotdf$sig != "", ], aes(color = .data$sig)) + 
        coord_fixed() + 
        scale_color_manual(values = genecolors[seq_len(2)], name = "") + 
        labs(x = substitute(Delta * "intron" ~ x, list(x = contrastName)),
             y = substitute(Delta * "exon" ~ x, list(x = contrastName))) + 
        theme_bw(base_size = 14) + 
        theme(panel.grid = element_blank(), 
              legend.position = "inside", 
              legend.position.inside = c(0.97, 0.03), 
              legend.justification = c(1, 0))
}
