#' Extract intron sequences
#'
#' @param gtf The path to a gtf file
#' @param genome A \code{DNAStringSet} object with the genome sequence
#' @param type Either 'collapse' or 'separate'
#' @param flanklength The length of the exonic flanking sequence
#'
#' @return A \code{DNAStringSet} object with intronic sequences
#' @export
#'
extractIntronSeqs <- function(gtf, genome, type = "collapse", flanklength = 90,
                              joinOverlappingIntrons = FALSE) {
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
    grl <- BiocGenerics::setdiff(range(grl), grl)
    
    ## Here, the order of the introns doesn't really matter, since they
    ## will not be joined together into a transcript
    
    ## Add flanking region
    grl <- grl + flanklength
    
    if (joinOverlappingIntrons) {
        ## If two (introns + flanklength) overlap, join them
        grl <- GenomicRanges::reduce(grl)
    }

    gr <- unlist(grl)
    
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

