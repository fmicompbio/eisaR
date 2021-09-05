## Genome
genome <- Biostrings::DNAStringSet(
    c(chr1 = "TTAACATTCGCTGGGGGAGATGACGAGACTAGCCGCCGCGTGGTCCTGCCGCATTATACGTGTTCAAGCGCCTACGTGGGTTGGGCAACCCGTGCCTATGGAGGCATGGACAAATTAGGTTCAACTTCAGCTACGTACGAGACCTAGAGGTAATAAGGGTATTTTACTCGGAGCATGTTTCAGTACGAACGTTAGATATC",
      chr2 = "CTATCGAAGTGGAATCTTGAAGAGCCCATCGGTTAAGGTCTCTCCAATGTCCAGCCTATTCTATGGCACGGCAGACCCGTTGTGCATCCACAGTGATAACTTACTTGGGCTCTTAATAGAGGAGTGTTGCCATTTTATCGGCTTGCACTCCAATTAGCACCAAGTGCCGTTATTGGGGTATTGCACTCATCAATAGCGTG")
)
genome_revcompl <- Biostrings::reverseComplement(genome)

## Transcript-to-gene mapping
t2g <- data.frame(
    transcript_id = c("tx1.1", "tx1.2", "tx1.3", "tx2.1-I", "tx2.2-U", "tx3.1", "tx3.2", "tx3.3"),
    gene_id = c("g1", "g1", "g1", "g2-I", "g2-I", "g3-U", "g3-U", "g3-U"),
    stringsAsFactors = FALSE
)

## Transcripts
txome <- Biostrings::DNAStringSet(
    c(`tx1.1` = paste0(substr(genome[["chr1"]], 116, 120),
                       substr(genome[["chr1"]], 141, 160)),
      `tx1.2` = paste0(substr(genome[["chr1"]], 61, 80),
                       substr(genome[["chr1"]], 101, 120),
                       substr(genome[["chr1"]], 141, 160)),
      `tx1.3` = paste0(substr(genome[["chr1"]], 61, 90)),
      `tx2.1-I` = paste0(substr(genome_revcompl[["chr1"]], 131, 150),
                         substr(genome_revcompl[["chr1"]], 171, 190)),
      `tx2.2-U` = paste0(substr(genome_revcompl[["chr1"]], 121, 140),
                         substr(genome_revcompl[["chr1"]], 156, 165),
                         substr(genome_revcompl[["chr1"]], 181, 190)),
      `tx3.1` = paste0(substr(genome[["chr2"]], 61, 90),
                       substr(genome[["chr2"]], 141, 160)),
      `tx3.2` = paste0(substr(genome[["chr2"]], 61, 80),
                       substr(genome[["chr2"]], 101, 120),
                       substr(genome[["chr2"]], 141, 160)),
      `tx3.3` = paste0(substr(genome[["chr2"]], 61, 90))
    )
)

## GTF
gtf <- GenomicRanges::GRanges(
    seqnames = rep(c("chr1", "chr2"), c(18, 10)),
    ranges = IRanges::IRanges(
        start = c(116, 141, 61, 101, 141, 61, 116, 61, 61, 61, 
                  51, 11, 61, 36, 11, 11, 11, 11, 
                  61, 141, 61, 101, 141, 61, 61, 61, 61, 61),
        end = c(120, 160, 80, 120, 160, 90, 160, 160, 90, 160, 
                70, 30, 80, 45, 20, 70, 80, 80, 
                90, 160, 80, 120, 160, 90, 160, 160, 90, 160)
    ),
    strand = rep(c("+", "-", "+"), c(10, 8, 10)),
    type = c(rep("exon", 6), rep("transcript", 3), "gene", rep("exon", 5), rep("transcript", 2), "gene",
             rep("exon", 6), rep("transcript", 3), "gene"),
    gene_id = rep(c("g1", "g2-I", "g3-U"), c(10, 8, 10)),
    transcript_id = c(rep("tx1.1", 2), rep("tx1.2", 3), "tx1.3", "tx1.1", "tx1.2", "tx1.3", NA, 
                      rep("tx2.1-I", 2), rep("tx2.2-U", 3), "tx2.1-I", "tx2.2-U", NA,
                      rep("tx3.1", 2), rep("tx3.2", 3), "tx3.3", "tx3.1", "tx3.2", "tx3.3", NA),
    exon_id = c("E1", "E2", "E3", "E4", "E5", "E6", NA, NA, NA, NA,
                "E7", "E8", "E9", "E10", "E11", NA, NA, NA,
                "E12", "E13", "E14", "E15", "E16", "E17", NA, NA, NA, NA)
)

rtracklayer::export(gtf, con = "small_example.gtf", format = "gtf")
Biostrings::writeXStringSet(txome, filepath = "small_example_txome.fa")
Biostrings::writeXStringSet(genome, filepath = "small_example_genome.fa")
write.table(t2g, file = "small_example_t2g.txt", row.names = FALSE,
            col.names = FALSE, quote = FALSE, sep = "\t")
