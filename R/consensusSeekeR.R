#' consensusSeekeR: Detection of consensus peak regions inside a group of
#' experiments using narrowPeak files
#'
#' This package compares narrowPeak data from multiple experiments to extract
#' common consensus peak regions. The size of the analyzed region is adjustable
#' as well as the number of experiments in which a peak must be detected
#' to mark a potential region as a consensus peak region.
#'
#' @docType package
#' @name consensusSeekeR-package
#' @aliases consensusSeekeR-package consensusSeekeR
#' @author  Astrid Louise Deschenes,
#'  Fabien Claude Lamaze,
#'  Pascal Belleau and
#'  Arnaud Droit
#'
#' Maintainer: Astrid Louise Deschenes <astrid-louise.deschenes@@crchudequebec.ulaval.ca>
#' @seealso
#'  \itemize{
#'    \item \code{\link{readNarrowPeakFile}} {for extracting regions and peaks
#'                  from a narrowPeak file.}
#'    \item \code{\link{findConsensusPeakRegions}} { for extracting regions
#'                  sharing the same features in more than one experiment. }
#'    }
#' @keywords package
NULL

#' Genomic regions with the greatest evidence of transcription factor binding
#' for the FOSL2 transcription factor (for demonstration purpose)
#'
#' Genomic regions representing the greatest evidence of enrichment for
#' the FOSL2 transcription factor (DCC accession: ENCFF002CFN)
#' for regions chr1:249120200-249250621 and chr10:1-370100
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' @name A549_FOSL2_01_NarrowPeaks_partial
#' @docType data
#' @format A \code{GRanges} containing one entry per genomic regions. Each row
#'  of \code{GRanges} has a name which represent the name of the experiment.
#' @references
#'  \itemize{
#'  \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia of
#'  DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#'  }
#' @seealso
#'  \itemize{
#'    \item \code{\link{A549_FOSL2_01_Peaks_partial}} { the associate
#'                  peaks dataset.}
#'    \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @usage data(A549_FOSL2_01_NarrowPeaks_partial)
#' @keywords datasets
#' @examples
#' ## Loading datasets
#' data(A549_FOSL2_01_NarrowPeaks_partial)
#' data(A549_FOXA1_01_NarrowPeaks_partial)
#' data(A549_FOSL2_01_Peaks_partial)
#' data(A549_FOXA1_01_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical name
#' names(A549_FOXA1_01_Peaks_partial) <- rep("FOXA1_01",
#'                              length(A549_FOXA1_01_Peaks_partial))
#' names(A549_FOXA1_01_NarrowPeaks_partial) <- rep("FOXA1_01",
#'                              length(A549_FOXA1_01_NarrowPeaks_partial))
#' names(A549_FOSL2_01_Peaks_partial) <-rep("FOSL2_01",
#'                              length(A549_FOSL2_01_Peaks_partial))
#' names(A549_FOSL2_01_NarrowPeaks_partial) <- rep("FOSL2_01",
#'                              length(A549_FOSL2_01_NarrowPeaks_partial))
#'
#' ## Calculating consensus regions for chromosome 10
#' chrList <- Seqinfo("chr10", 135534747, NA)
#' findConsensusPeakRegions(
#'      narrowPeaks = c(A549_FOXA1_01_NarrowPeaks_partial, A549_FOSL2_01_NarrowPeaks_partial),
#'      peaks = c(A549_FOXA1_01_Peaks_partial, A549_FOSL2_01_Peaks_partial),
#'      chrInfo = chrList, extendingSize = 100, includeAllPeakRegion = FALSE,
#'      minNbrExp = 2, nbrThreads = 1)
#'
NULL

#' TODO Data Set
#'
#' TODO
#'
#' @name A549_FOSL2_01_Peaks_partial
#' @docType data
#' @format A \code{GRanges}
#' @usage data(A549_FOSL2_01_Peaks_partial)
#' @references
#'  \itemize{
#'  \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia of
#'  DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#'  }
#' @keywords datasets
NULL

#' TODO Data Set
#'
#' TODO
#'
#' @name A549_FOXA1_01_NarrowPeaks_partial
#' @docType data
#' @format A \code{GRanges}
#' @usage data(A549_FOXA1_01_NarrowPeaks_partial)
#' @references
#'  \itemize{
#'  \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia of
#'  DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#'  }
#' @keywords datasets
NULL

#' TODO Data Set
#'
#' @name A549_FOXA1_01_Peaks_partial
#' @docType data
#' @format A \code{GRanges}
#' @usage data(A549_FOXA1_01_Peaks_partial)
#' @references
#'  \itemize{
#'  \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia of
#'  DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#'  }
#' @keywords datasets
NULL