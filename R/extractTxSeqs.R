#' Extract (spliced or unspliced) transcript sequences
#'
#' @param gtf The path to a gtf file
#' @param genome A \code{DNAStringSet} object with the genome sequence
#' @param type Either 'spliced' or 'unspliced'
#'
#' @return A \code{DNAStringSet} object with intronic sequences
#' @export
#'
#' @importFrom BSgenome getSeq
#' @importFrom GenomicFeatures extractTranscriptSeqs makeTxDbFromGFF exonsBy
#'
extractTxSeqs <- function(gtf, genome, type = "spliced") {
    ## TODO: Check that gtf is path to existing file, genome is DNAStringSet, type is either 'spliced' or 'unspliced'
    
    ## Construct TxDb from gtf file
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")
    
    ## Group exons by transcript. At this stage, exons are ordered
    ## by exon_rank (i.e., from 5' to 3' in the transcript)
    grl <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
    
    ## TODO: Check that exons are ordered by exon_rank
    
    ## Extract transcript sequences
    if (type == "spliced") {
        txout <- GenomicFeatures::extractTranscriptSeqs(x = genome, transcripts = grl)
    } else if (type == "unspliced") {
        txout <- GenomicFeatures::extractTranscriptSeqs(x = genome, transcripts = range(grl))
    } else {
        stop("Unknown 'type' argument")
    }
    
    return(txout)
}
