# eisaR 1.25.1

* Rewrite plotEISA to use ggplot2 (and return a ggplot object), set the aspect 
  ratio of the plot to 1 and add an identity line. Additional arguments 
  provided via ..., previously passed to plot(), are ignored.

# eisaR 1.21.1

* Various code improvements detected by lintr

# eisaR 1.15.1

* Add legacyQLF argument to runEISA (will be passed to edgeR::glmQLFit)

# eisaR 1.5.2

* Add options to collapse introns by gene and restrict introns to feature 
  ranges in getFeatureRanges.

# eisaR 1.0

* Initial version  
* R implementation of Exon-Intron Split Analysis (EISA, doi: 10.1038/nbt.3269)  
* functionality for extracting spliced and unspliced transcript sequences,
  as well as intron sequences. e.g. for use in RNA velocity analysis
