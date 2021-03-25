#' @title Get exonic/gene body regions from a transcript database.
#'
#' @description From a transcript database package (\code{\link[GenomicFeatures:TxDb-class]{TxDb}}),
#'   extract exonic and gene body ranges for use with EISA. These regions can
#'   be used to quantify RNA-seq alignments in exons and gene bodies, respectively.
#'   Intronic counts can then be obtained from the difference between gene bodies
#'   and exonic region counts.
#'
#' @author Michael Stadler
#'
#' @param txdb a \code{TxDb} or an \code{EnsDb} object with the transcript annotations.
#' @param exonExt \code{numeric} (default = 10L). Exonic ranges will be extended
#'   on either side by this many nucleotides, in order to avoid "bleed-over" of
#'   exonic alignments into adjacent intronic regions.
#' @param strandedData \code{logical(1)}. If \code{TRUE}, the RNA-seq data is
#'   assumed to be strand-specific, and therefore only overlapping genes that
#'   are on the same strand will be filtered out. If \code{FALSE}, also genes
#'   overlapping on opposite strands will be filtered out.
#'
#' @details The exonic regions are generated as follows:
#'   \enumerate{
#'     \item extract exons by gene from the \code{txdb}
#'     \item extend each exon by \code{exonExt}
#'     \item combine overlapping exons within each gene
#'     \item create gene body ranges from the most extreme exonic coordinates
#'     \item filter out genes that have only a single exon (no intron), have exons
#'           on more than a single chromosome or on both strands, or that
#'           overlap other genes
#'     }
#' @return a \code{list} with elements "exons" and "genebodies", containing
#'   named \code{GenomicRanges} objects with ranges for exons and gene bodies,
#'   respectively.
#'
#' @seealso \code{\link[GenomicFeatures:TxDb-class]{TxDb}} for details on \code{TxDb} objects
#'   and how to create them, e.g. from \code{.gtf} files.
#'
#' @examples
#' if (requireNamespace("AnnotationDbi", quietly = TRUE)) {
#'     txdb <- AnnotationDbi::loadDb(system.file("extdata", "hg19sub.sqlite", package = "eisaR"))
#'     regL <- getRegionsFromTxDb(txdb)
#'     lengths(regL)
#' }
#'
#' @import GenomicRanges
#' @importFrom S4Vectors elementNROWS runLength
#' @importFrom IRanges IRanges
#'
#' @export
getRegionsFromTxDb <- function(txdb, exonExt = 10L, strandedData = TRUE) {
    # check arguments
    stopifnot(inherits(txdb, "TxDb") || inherits(txdb, "EnsDb"))
    stopifnot(is.numeric(exonExt) && length(exonExt) == 1L)
    stopifnot(is.logical(strandedData) && length(strandedData) == 1L)
    for (pkg in c("GenomicFeatures", "AnnotationDbi")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("getFeatureRanges() requires installing the Bioconductor package '",
                 pkg, "' using BiocManager::install(\"", pkg, "\")")
        }
    }
    if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
        stop("getRegionsFromTxDb() requires installing the Bioconductor package 'GenomicFeatures'",
             " using BiocManager::install(\"GenomicFeatures\")")
    }
    
    # get exons and extend by exonExt on both sides
    message("extracting exon coordinates")
    exL <- GenomicFeatures::exonsBy(txdb, by = "gene")
    suppressWarnings(exL <- exL + exonExt) # may contain out of chromosome regions
    # remark: should GenomicRanges::trim now, but there is not trim,GRangesList method
    #         will suppressWarnings below until I can use trim,GRanges
    exL <- GenomicRanges::reduce(exL) # fuse

    # identify genes with single exon
    nEx <- S4Vectors::elementNROWS(exL)
    severalExons <- nEx > 1

    # identify genes with exons on multiple chromosomes
    suppressWarnings(nChr <- S4Vectors::elementNROWS(S4Vectors::runLength(GenomicRanges::seqnames(exL))))
    singleChr <- nChr == 1

    ## identify genes with exons on multiple strands
    suppressWarnings(nStr <- S4Vectors::elementNROWS(S4Vectors::runLength(GenomicRanges::strand(exL))))
    singleStrand <- nStr == 1

    # get gene body
    ex <- suppressWarnings(GenomicRanges::trim(unlist(exL)))
    suppressWarnings(gbody <- unlist(range(exL))[unique(names(ex))])

    # trim out-of chromosome regions
    gbody <- GenomicRanges::trim(gbody)

    ## identify overlapping genes
    nov <- GenomicRanges::countOverlaps(gbody, gbody, ignore.strand = !strandedData)
    noOverlaps <- nov == 1 # only self-overlaps

    ## filter regions
    sel <- severalExons & singleChr & singleStrand & noOverlaps
    selex <- names(ex) %in% names(gbody[sel])
    message("total number of genes/exons: ", length(sel), "/", length(ex))
    message("removing overlapping/single-exon/ambiguous genes (", length(sel) - sum(sel), ")")
    message("creating filtered regions for ", sum(sel), " genes (", round(mean(sel)*100, 1),
            "%) with ", sum(selex), " exons (", round(mean(selex) * 100, 1), "%)")
    gbody <- gbody[sel]
    ex <- ex[selex]

    ## return results
    return(list(exons = ex, genebodies = gbody))
}

