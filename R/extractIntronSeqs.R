#' Extract intron sequences
#'
#' @param gtf The path to a gtf file
#' @param genome A \code{DNAStringSet} object with the genome sequence
#' @param type Either 'collapse' or 'separate'
#' @param flanklength The length of the exonic/genomic flanking sequence that 
#'   is added to each intron
#' @param joinOverlappingIntrons Whether to join introns of the same gene if
#'   they overlap after expanding with the flanking sequences
#' @param verbose Logical, whether or not to print messages
#'
#' @return A \code{DNAStringSet} object with intronic sequences
#' @export
#' @author Charlotte Soneson
#' 
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' @importFrom GenomicRanges reduce
#' @importFrom BiocGenerics setdiff unlist
#' @importFrom BSgenome getSeq
#'
extractIntronSeqs <- function(gtf, genome, type = c("collapse", "separate"), 
                              flanklength = 90L, joinOverlappingIntrons = FALSE,
                              verbose = TRUE) {
    type <- match.arg(type)
    
    if (!is.character(gtf) || length(gtf) != 1 || !file.exists(gtf)) {
        stop("'gtf' must be a character scalar providing ", 
             "the path to an existing file.")
    }
    if (!is(genome, "DNAStringSet")) {
        stop("'genome' must be a DNAStringSet")
    }
    if (!is.numeric(flanklength) || length(flanklength) != 1 ||
        flanklength <= 0) {
        stop("'flanklength' must be a numeric non-negative scalar")
    }
    flanklength <- as.integer(flanklength)
    if (!is.logical(joinOverlappingIntrons) || length(joinOverlappingIntrons) != 1) {
        stop("'joinOverlappingIntrons' must be a logical scalar")
    }
    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be a logical scalar")
    }
    
    if (verbose) {
        message("Extracting intron sequences using the ", type, " approach")
    }
    
    ## Construct TxDb from gtf file
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")
    
    if (type == "separate") {
        ## Group exons by transcript
        grl <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
    } else if (type == "collapse") {
        ## Group exons by gene
        grl <- GenomicFeatures::exonsBy(txdb, by = "gene")
        
        ## Collapse the exons of each gene
        grl <- GenomicRanges::reduce(grl)
    } else {
        stop("Unknown 'type' argument")
    }
    
    ## Get introns as the set difference between the range and the exons,
    ## for each transcript
    ## Here, the order of the introns doesn't really matter, since they
    ## will not be joined together into a transcript
    grl <- BiocGenerics::setdiff(range(grl), grl)
    
    ## Add flanking region
    grl <- grl + flanklength
    
    if (joinOverlappingIntrons) {
        ## If two (introns + flanklength) overlap, join them
        grl <- GenomicRanges::reduce(grl)
    }

    gr <- BiocGenerics::unlist(grl)
    
    ## Add -I{X} to names
    names(gr) <- gsub("\\-I\\.", "-I", make.unique(paste0(names(gr), "-I")))
    
    ## Get sequence
    gs <- BSgenome::getSeq(x = genome, names = gr)
    
    ## Manually set names of extracted sequences
    if (any(width(gs) != width(gr))) {
        stop("Something went wrong in the sequence extraction - ", 
             "the widths of the extracted sequences don't agree ", 
             "with the widths of the features")
    }
    if (any(names(gs) != names(gr))) {
        warning("Setting names of extracted sequences manually")
        names(gs) <- names(gr)
    }

    return(gs)
}

