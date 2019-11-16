context("Reference generation")

test_that("reference generation works", {
    ## Generate simple example
    genome <- Biostrings::DNAStringSet(c(chr1 = "CAACACCGACTGGATTTGAGATGACCGTACCGATCCAAATCTCAGTGGGGCAACATCAGCCTTATAACGCTCGACTTTGACAATTATTCTGTCCGCTCTA",
                                         chr2 = "CAACACCGACTGGATTTGAGATGACCGTACCGATCCAAATCTCAGTGGGGCAACATCAGCCTTATAACGCTCGACTTTGACAATTATTCTGTCCGCTCTA"))
    genome_rev <- Biostrings::reverseComplement(genome)
    txome <- Biostrings::DNAStringSet(
        c(g1.tx1 = paste0(substr(genome[["chr1"]], 1, 20),
                          substr(genome[["chr1"]], 81, 100)),
          g1.tx2 = paste0(substr(genome[["chr1"]], 1, 30),
                          substr(genome[["chr1"]], 81, 100)),
          g1.tx3 = paste0(substr(genome[["chr1"]], 1, 20),
                          substr(genome[["chr1"]], 41, 60),
                          substr(genome[["chr1"]], 81, 100)),
          g2.tx1 = paste0(substr(genome_rev[["chr2"]], 1, 20),
                          substr(genome_rev[["chr2"]], 81, 100)),
          g2.tx2 = paste0(substr(genome_rev[["chr2"]], 1, 20),
                          substr(genome_rev[["chr2"]], 71, 100)),
          g2.tx3 = paste0(substr(genome_rev[["chr2"]], 1, 20),
                          substr(genome_rev[["chr2"]], 41, 60),
                          substr(genome_rev[["chr2"]], 81, 100)))
    )
    gtf <- GenomicRanges::GRanges(
        seqnames = rep(c("chr1", "chr2"), each = 11),
        ranges = IRanges::IRanges(start = c(1, 81, 1, 81, 1, 41, 81, 1, 1, 1, 1,
                                            81, 1, 1, 81, 81, 41, 1, 1, 1, 1, 1),
                                  end = c(20, 100, 30, 100, 20, 60, 100, 100, 100, 100, 100,
                                          100, 20, 30, 100, 100, 60, 20, 100, 100, 100, 100)),
        strand = rep(c("+", "-"), each = 11),
        type = rep(c(rep("exon", 7), rep("transcript", 3), "gene"), 2),
        gene_id = rep(c("g1", "g2"), each = 11),
        transcript_id = c("g1.tx1", "g1.tx1", "g1.tx2", "g1.tx2", "g1.tx3", "g1.tx3",
                          "g1.tx3", "g1.tx1", "g1.tx2", "g1.tx3", NA,
                          "g2.tx1", "g2.tx1", "g2.tx2", "g2.tx2", "g2.tx3", "g2.tx3",
                          "g2.tx3", "g2.tx1", "g2.tx2", "g2.tx3", NA),
        exon_id = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", NA, NA, NA, NA,
                    "E8", "E9", "E10", "E11", "E12", "E13", "E14", NA, NA, NA, NA)
    )
    f <- tempfile()
    rtracklayer::export(gtf, con = f, format = "gtf")
    
    txs <- extractTxSeqs(gtf = f, genome = genome, type = "spliced")
    expect_equal(as.character(txs), as.character(txome))
    
    txu <- extractTxSeqs(gtf = f, genome = genome, type = "unspliced")
    expect_equal(as.character(txu[["g1.tx1"]]), as.character(genome[["chr1"]]))
    expect_equal(as.character(txu[["g1.tx2"]]), as.character(genome[["chr1"]]))
    expect_equal(as.character(txu[["g1.tx2"]]), as.character(genome[["chr1"]]))
    expect_equal(as.character(txu[["g2.tx1"]]), as.character(genome_rev[["chr2"]]))
    expect_equal(as.character(txu[["g2.tx2"]]), as.character(genome_rev[["chr2"]]))
    expect_equal(as.character(txu[["g2.tx3"]]), as.character(genome_rev[["chr2"]]))
    
    inc <- extractIntronSeqs(gtf = f, genome = genome, type = "collapse",
                             flanklength = 3)
    expect_equal(as.character(inc[["g1-I"]]),
                 as.character(substr(genome["chr1"], 28, 43)))
    expect_equal(as.character(inc[["g1-I1"]]),
                 as.character(substr(genome["chr1"], 58, 83)))
    expect_equal(as.character(inc[["g2-I"]]),
                 as.character(substr(genome_rev["chr2"], 58, 73)))
    expect_equal(as.character(inc[["g2-I1"]]),
                 as.character(substr(genome_rev["chr2"], 18, 43)))
    
    ins <- extractIntronSeqs(gtf = f, genome = genome, type = "separate",
                             flanklength = 3)
    expect_equal(as.character(ins[["g1.tx1-I"]]),
                 as.character(substr(genome["chr1"], 18, 83)))
    expect_equal(as.character(ins[["g1.tx2-I"]]),
                 as.character(substr(genome["chr1"], 28, 83)))
    expect_equal(as.character(ins[["g1.tx3-I"]]),
                 as.character(substr(genome["chr1"], 18, 43)))
    expect_equal(as.character(ins[["g1.tx3-I1"]]),
                 as.character(substr(genome["chr1"], 58, 83)))
    expect_equal(as.character(ins[["g2.tx1-I"]]),
                 as.character(substr(genome_rev["chr2"], 18, 83)))
    expect_equal(as.character(ins[["g2.tx2-I"]]),
                 as.character(substr(genome_rev["chr2"], 18, 73)))
    expect_equal(as.character(ins[["g2.tx3-I"]]),
                 as.character(substr(genome_rev["chr2"], 58, 83)))
    expect_equal(as.character(ins[["g2.tx3-I1"]]),
                 as.character(substr(genome_rev["chr2"], 18, 43)))
    
})
