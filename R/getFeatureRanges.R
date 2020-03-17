#' Generate a GRangesList object with feature ranges
#' 
#' Generate a GRangesList object with genomic ranges for (any combination of) 
#' spliced transcripts, unspliced transcripts and introns.
#' 
#' @param gtf Path to gtf file.
#' @param featureType Character vector indicating the type(s) of features to
#'   extract, any subset of \code{c("spliced", "intron", "unspliced")}.
#' @param intronType Character vector indicating how to define the introns 
#'   (only used if "intron" is part of \code{featureType}). Has to be either 
#'   "separate" (introns are defined for each transcript separately) or 
#'   "collapse" (transcripts of the same gene are first collapsed before 
#'   introns are defined as any non-exonic part of the gene locus). 
#' @param flankLength Integer scalar indicating the length of the flanking 
#'   sequence added to each side of each extracted intron (only used if 
#'   "intron" is included in \code{featureType}).
#' @param joinOverlappingIntrons Logical scalar indicating whether two introns 
#'   that overlap (after adding the flanking sequences) should be joined into 
#'   one feature.
#' @param verbose Logical scalar, whether to print out progress messages.
#' 
#' @author Charlotte Soneson
#' 
#' @export
#' 
#' @return Returns a \code{GRangesList} object where each element represents 
#'   one extracted feature. The metadata of this object contains two 
#'   \code{data.frame}s mapping corresponding identifiers between the 
#'   different feature types, as well as a list of all features for each type.
#' 
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' @importFrom GenomicRanges GRangesList reduce
#' @importFrom AnnotationDbi select
#' @importFrom BiocGenerics unlist relist setdiff
#' @importFrom S4Vectors mcols metadata
#' 
#' @examples
#'   ## Get feature ranges
#'   grl <- getFeatureRanges(
#'     gtf = system.file("extdata/small_example.gtf", package = "eisaR"),
#'     featureType = c("spliced", "intron"),
#'     intronType = "separate",
#'     flankLength = 5L,
#'     joinOverlappingIntrons = FALSE,
#'     verbose = TRUE
#'   )
#'   
#'   ## GRangesList
#'   grl
#'   
#'   ## Corresponding transcript/gene IDs
#'   S4Vectors::metadata(grl)$corrtx
#'   S4Vectors::metadata(grl)$corrgene
#'   
#'   ## List of features of different types
#'   S4Vectors::metadata(grl)$featurelist
#'   
#'   ## Get feature sequences
#'   if (requireNamespace("BSgenome", quietly = TRUE)) {
#'     library(BSgenome)
#'     genome <- Biostrings::readDNAStringSet(
#'       system.file("extdata/small_example_genome.fa", package = "eisaR"))
#'     seqs <- GenomicFeatures::extractTranscriptSeqs(x = genome, 
#'                                                    transcripts = grl)
#'     seqs
#'   }
#'  
getFeatureRanges <- function(
    gtf, 
    featureType = c("spliced", "intron"), 
    intronType = "separate", 
    flankLength = 90L, 
    joinOverlappingIntrons = FALSE, 
    verbose = TRUE) {
    
    ## --------------------------------------------------------------------- ##
    ## Pre-flight checks
    ## --------------------------------------------------------------------- ##
    if (length(gtf) != 1 || !is.character(gtf) || !file.exists(gtf)) {
        stop("'gtf' must be a character scalar providing ", 
             "the path to an existing file.")
    }
    if (!all(is.character(featureType)) || 
        !all(featureType %in% c("spliced", "unspliced", "intron"))) {
        stop("'featureType' must be a subset of c('spliced', 'unspliced', 'intron')")
    }
    if ("intron" %in% featureType) {
        if (length(intronType) != 1 || !is.character(intronType) || 
            !(intronType %in% c("separate", "collapse"))) {
            stop("'intronType' must be either 'separate' or 'collapse'")
        }
        if (!is.numeric(flankLength) || length(flankLength) != 1 ||
            flankLength < 0) {
            stop("'flankLength' must be a numeric non-negative scalar")
        }
        flankLength <- as.integer(flankLength)
        if (!is.logical(joinOverlappingIntrons) || 
            length(joinOverlappingIntrons) != 1) {
            stop("'joinOverlappingIntrons' must be a logical scalar")
        }
    }
    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("'verbose' must be a logical scalar")
    }
    
    ## --------------------------------------------------------------------- ##
    ## Define suffixes for unspliced transcripts and introns
    ## --------------------------------------------------------------------- ##
    suffixes <- c(intron = "-I", unspliced = "-U")
    
    ## --------------------------------------------------------------------- ##
    ## Construct TxDb from gtf file
    ## --------------------------------------------------------------------- ##
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")
    
    ## Initialize GRangesList
    grlfull <- GenomicRanges::GRangesList()
    featurelist <- list()
    
    ## --------------------------------------------------------------------- ##
    ## Group exons by transcript/gene
    ## --------------------------------------------------------------------- ##
    ## Group exons by transcript. When using exonsBy with by = "tx", 
    ## the returned exons are ordered by ascending rank for each transcript, 
    ## that is, by their position in the transcript. 
    ebt <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
    t2g <- AnnotationDbi::select(txdb, keys = names(ebt), 
                                 keytype = "TXNAME", columns = "GENEID")
    e2 <- BiocGenerics::unlist(ebt)
    e2$transcript_id <- names(e2)
    e2$gene_id = t2g$GENEID[match(e2$transcript_id, t2g$TXNAME)]
    e2$exon_id <- e2$exon_name
    e2$exon_name <- NULL
    e2$type <- "exon"
    names(e2) <- NULL
    mcols(e2) <- S4Vectors::mcols(e2)[, c("exon_id", "exon_rank", 
                                          "transcript_id", "gene_id", "type")]
    ebt <- BiocGenerics::relist(e2, ebt)
    
    ## Group exons by gene
    ebg <- GenomicFeatures::exonsBy(txdb, by = "gene")
    
    corrtx <- data.frame(spliced = unique(names(ebt)),
                         stringsAsFactors = FALSE)
    corrtx$unspliced <- paste0(corrtx$spliced, suffixes["unspliced"])
    corrgene <- data.frame(spliced = unique(unlist(ebt)$gene_id),
                           stringsAsFactors = FALSE)
    corrgene$unspliced <- paste0(corrgene$spliced, suffixes["unspliced"])
    corrgene$intron <- paste0(corrgene$spliced, suffixes["intron"])
    
    ## --------------------------------------------------------------------- ##
    ## Spliced transcripts
    ## --------------------------------------------------------------------- ##
    if ("spliced" %in% featureType) {
        if (verbose) {
            message("Extracting spliced transcript features")
        }
        
        featurelist$spliced <- names(ebt)
        
        ## Here, it's important that for each transcript, the exons must be ordered 
        ## by ascending rank, that is, by ascending position in the transcript.
        grlfull <- c(grlfull, ebt)
    }
    
    ## --------------------------------------------------------------------- ##
    ## Unspliced transcripts
    ## --------------------------------------------------------------------- ##
    if ("unspliced" %in% featureType) {
        if (verbose) {
            message("Extracting unspliced transcript features")
        }
        
        ebtr <- range(ebt)
        e2 <- BiocGenerics::unlist(ebtr)
        e2$exon_rank <- 1L
        e2$transcript_id <- names(e2)
        e2$gene_id = t2g$GENEID[match(e2$transcript_id, t2g$TXNAME)]
        e2$type <- "exon"
        e2$transcript_id <- paste0(e2$transcript_id, suffixes["unspliced"])
        e2$gene_id <- paste0(e2$gene_id, suffixes["unspliced"])
        e2$exon_id <- e2$transcript_id
        names(e2) <- NULL
        mcols(e2) <- S4Vectors::mcols(e2)[, c("exon_id", "exon_rank", 
                                              "transcript_id", "gene_id", "type")]
        ebtr <- BiocGenerics::relist(e2, ebtr)
        names(ebtr) <- paste0(names(ebtr), suffixes["unspliced"])
        
        featurelist$unspliced <- names(ebtr
                                       )
        grlfull <- c(grlfull, ebtr)
    }
    
    ## --------------------------------------------------------------------- ##
    ## Introns
    ## --------------------------------------------------------------------- ##
    if ("intron" %in% featureType) {
        if (verbose) {
            message("Extracting introns using the ", intronType, " approach")
        }
        
        if (intronType == "separate") {
            ## Group exons by transcript
            grl <- ebt
        } else if (intronType == "collapse") {
            ## Group exons by gene
            grl <- ebg
            
            ## Collapse the exons of each gene
            grl <- GenomicRanges::reduce(grl)
        }
        
        ## Get introns as the set difference between the range and the exons,
        ## for each transcript
        ## Here, the order of the introns doesn't really matter, since they
        ## will not be joined together into a transcript
        grl <- BiocGenerics::setdiff(range(grl), grl)
        
        ## Remove empty entries
        grl <- grl[vapply(grl, length, 0L) > 0]
        
        ## Add flanking region
        grl <- grl + flankLength
        
        ## Make sure that introns don't go to negative coordinates
        grl <- GenomicRanges::restrict(grl, start = 1L)
        
        if (joinOverlappingIntrons) {
            ## If two (introns + flankLength) overlap, join them
            grl <- GenomicRanges::reduce(grl)
        }

        gr <- BiocGenerics::unlist(grl)
        gr$exon_rank <- 1L
        gr$transcript_id <- names(gr)
        if (intronType == "separate") {
            gr$gene_id <- t2g$GENEID[match(gr$transcript_id, t2g$TXNAME)]
        } else {
            gr$gene_id <- gr$transcript_id
        }
        gr$type <- "exon"
        gr$transcript_id <- gsub(paste0(suffixes["intron"], "."), suffixes["intron"], 
                                 make.unique(paste0(gr$transcript_id, suffixes["intron"])),
                                 fixed = TRUE)
        gr$gene_id <- paste0(gr$gene_id, suffixes["intron"])
        gr$exon_id <- gr$transcript_id
        names(gr) <- NULL
        mcols(gr) <- S4Vectors::mcols(gr)[, c("exon_id", "exon_rank", 
                                              "transcript_id", "gene_id", "type")]
        grl <- BiocGenerics::relist(gr, lapply(
            structure(seq_along(gr), names = gr$transcript_id), function(i) i))
        
        featurelist$intron <- names(grl)

        grlfull <- c(grlfull, grl)
    }
    
    S4Vectors::metadata(grlfull)$corrtx <- 
        corrtx[, colnames(corrtx) %in% featureType, drop = FALSE]
    S4Vectors::metadata(grlfull)$corrgene <- 
        corrgene[, colnames(corrgene) %in% featureType, drop = FALSE]
    S4Vectors::metadata(grlfull)$featurelist <- featurelist
    
    grlfull
}





