      #' consensusSeekeR: Detection of consensus peak regions inside a group of
#' experiments using narrowPeak files
#'
#' This package compares narrowPeak data from multiple experiments to extract
#' common consensus peak regions. The size of the analyzed region is adjustable
#' as well as the number of experiments in which a peak must be detected
#' to mark a potential region as a consensus peak region.
#'
#' @docType package
#'
#' @name consensusSeekeR-package
#'
#' @aliases consensusSeekeR-package consensusSeekeR
#'
#' @author  Astrid Louise Deschenes,
#' Fabien Claude Lamaze,
#' Pascal Belleau and
#' Arnaud Droit
#'
#' Maintainer:
#' Astrid Louise Deschenes <astrid-louise.deschenes@@crchudequebec.ulaval.ca>
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{readNarrowPeakFile}} {for extracting regions and peaks
#' from a narrowPeak file.}
#'     \item \code{\link{findConsensusPeakRegions}} { for extracting regions
#' sharing the same features in more than one experiment. }
#' }
#'
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
#'
#' @docType data
#'
#' @aliases A549_FOSL2_01_NarrowPeaks_partial
#'
#' @format A \code{GRanges} containing one entry per genomic regions. Each row
#' of \code{GRanges} has a name which represent the name of the experiment.
#'
#' @source The Encyclopedia of DNA Elements (ENCODE) (DCC accession:
#' ENCFF002CFN)
#'
#' @references
#' \itemize{
#'      \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_FOSL2_01_Peaks_partial}} { the associate
#' sites dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#' sharing the same features in more than one experiment.}
#' }
#'
#' @usage data(A549_FOSL2_01_NarrowPeaks_partial)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_FOSL2_01_NarrowPeaks_partial)
#' data(A549_FOSL2_01_Peaks_partial)
#' data(A549_FOXA1_01_NarrowPeaks_partial)
#' data(A549_FOXA1_01_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_FOXA1_01_Peaks_partial) <- rep("FOXA1_01",
#'                             length(A549_FOXA1_01_Peaks_partial))
#' names(A549_FOXA1_01_NarrowPeaks_partial) <- rep("FOXA1_01",
#'                             length(A549_FOXA1_01_NarrowPeaks_partial))
#' names(A549_FOSL2_01_Peaks_partial) <-rep("FOSL2_01",
#'                             length(A549_FOSL2_01_Peaks_partial))
#' names(A549_FOSL2_01_NarrowPeaks_partial) <- rep("FOSL2_01",
#'                             length(A549_FOSL2_01_NarrowPeaks_partial))
#'
#' ## Calculating consensus regions for chromosome 10 only
#' ## with a default region size of 200 bp (2 * extendingSize)
#' ## which is not extended to include all genomic regions.
#' ## A peak from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo("chr10", 135534747, NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_FOXA1_01_NarrowPeaks_partial,
#'                         A549_FOSL2_01_NarrowPeaks_partial),
#'     peaks = c(A549_FOXA1_01_Peaks_partial,
#'                         A549_FOSL2_01_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 100,
#'     expandToFitPeakRegion = FALSE,
#'     shrinkToFitPeakRegion = TRUE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL

#' Sites with the greatest evidence of transcription factor binding
#' for the FOSL2 transcription factor (for demonstration purpose)
#'
#' Sites representing the greatest evidence of enrichment for
#' the FOSL2 transcription factor (DCC accession: ENCFF002CFN)
#' for regions chr1:249120200-249250621 and chr10:1-370100
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' @name A549_FOSL2_01_Peaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per site. Each row
#' of \code{GRanges} has the same row name which represent the name
#' of the experiment.
#'
#' @usage data(A549_FOSL2_01_Peaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_FOSL2_01_NarrowPeaks_partial}} { the associate
#' genomic regions dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#' sharing the same features in more than one experiment.}
#' }
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_FOSL2_01_NarrowPeaks_partial)
#' data(A549_FOXA1_01_NarrowPeaks_partial)
#' data(A549_FOSL2_01_Peaks_partial)
#' data(A549_FOXA1_01_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_FOXA1_01_Peaks_partial) <- rep("FOXA1_01",
#'                             length(A549_FOXA1_01_Peaks_partial))
#' names(A549_FOXA1_01_NarrowPeaks_partial) <- rep("FOXA1_01",
#'                             length(A549_FOXA1_01_NarrowPeaks_partial))
#' names(A549_FOSL2_01_Peaks_partial) <-rep("FOSL2_01",
#'                             length(A549_FOSL2_01_Peaks_partial))
#' names(A549_FOSL2_01_NarrowPeaks_partial) <- rep("FOSL2_01",
#'                             length(A549_FOSL2_01_NarrowPeaks_partial))
#'
#' ## Calculating consensus regions for chromosome 1 only
#' ## with a default region size of 400 bp (2 * extendingSize)
#' ## which is extended to include all genomic regions of the
#' ## closest peak (for each experiment).
#' ## A peak from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo("chr1", 249250621, NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_FOXA1_01_NarrowPeaks_partial,
#'                         A549_FOSL2_01_NarrowPeaks_partial),
#'     peaks = c(A549_FOXA1_01_Peaks_partial,
#'                         A549_FOSL2_01_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 200,
#'     expandToFitPeakRegion = TRUE,
#'     shrinkToFitPeakRegion = FALSE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL

#' Genomic regions with the greatest evidence of transcription factor binding
#' for the FOXA1 transcription factor (for demonstration purpose)
#'
#' Genomic regions representing the greatest evidence of enrichment for
#' the FOXA1 transcription factor (DCC accession: TODO)
#' for regions chr1:249120200-249250621 and chr10:1-370100
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' @name A549_FOXA1_01_NarrowPeaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per genomic regions. Each row
#' of \code{GRanges} has a name which represent the name of the experiment.
#'
#' @usage data(A549_FOXA1_01_NarrowPeaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_FOXA1_01_Peaks_partial}} { the associate
#' sites dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#' sharing the same features in more than one experiment.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_FOSL2_01_NarrowPeaks_partial)
#' data(A549_FOSL2_01_Peaks_partial)
#' data(A549_FOXA1_01_NarrowPeaks_partial)
#' data(A549_FOXA1_01_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_FOXA1_01_Peaks_partial) <- rep("FOXA1_01",
#'                              length(A549_FOXA1_01_Peaks_partial))
#' names(A549_FOXA1_01_NarrowPeaks_partial) <- rep("FOXA1_01",
#'                              length(A549_FOXA1_01_NarrowPeaks_partial))
#' names(A549_FOSL2_01_Peaks_partial) <-rep("FOSL2_01",
#'                              length(A549_FOSL2_01_Peaks_partial))
#' names(A549_FOSL2_01_NarrowPeaks_partial) <- rep("FOSL2_01",
#'                              length(A549_FOSL2_01_NarrowPeaks_partial))
#'
#' ## Calculating consensus regions for both chromosomes 1 and 10
#' ## with a default region size of 300 bp (2 * extendingSize)
#' ## which is not extended to include all genomic regions.
#' ## A peak from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo(c("chr1", "chr10"), c(249250621, 135534747), NA)
#' findConsensusPeakRegions(
#'      narrowPeaks = c(A549_FOXA1_01_NarrowPeaks_partial,
#'                         A549_FOSL2_01_NarrowPeaks_partial),
#'      peaks = c(A549_FOXA1_01_Peaks_partial,
#'                         A549_FOSL2_01_Peaks_partial),
#'      chrInfo = chrList,
#'      extendingSize = 150,
#'      expandToFitPeakRegion = FALSE,
#'      minNbrExp = 2,
#'      nbrThreads = 1)
#'
NULL

#' Sites with the greatest evidence of transcription factor binding
#' for the FOXA1 transcription factor (for demonstration purpose)
#'
#' Sites representing the greatest evidence of enrichment for
#' the FOXA1 transcription factor (DCC accession: TODO)
#' for regions chr1:249120200-249250621 and chr10:1-370100
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' @name A549_FOXA1_01_Peaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per site
#' . Each row of \code{GRanges} has a name which represent the name of
#' the experiment.
#'
#' @usage data(A549_FOXA1_01_Peaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_FOXA1_01_NarrowPeaks_partial}} { the associate
#'                  genomic regions dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_FOSL2_01_NarrowPeaks_partial)
#' data(A549_FOSL2_01_Peaks_partial)
#' data(A549_FOXA1_01_NarrowPeaks_partial)
#' data(A549_FOXA1_01_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_FOXA1_01_Peaks_partial) <- rep("FOXA1_01",
#'                               length(A549_FOXA1_01_Peaks_partial))
#' names(A549_FOXA1_01_NarrowPeaks_partial) <- rep("FOXA1_01",
#'                               length(A549_FOXA1_01_NarrowPeaks_partial))
#' names(A549_FOSL2_01_Peaks_partial) <-rep("FOSL2_01",
#'                               length(A549_FOSL2_01_Peaks_partial))
#' names(A549_FOSL2_01_NarrowPeaks_partial) <- rep("FOSL2_01",
#'                               length(A549_FOSL2_01_NarrowPeaks_partial))
#'
#' ## Calculating consensus regions for both chromosomes 1 and 10
#' ## with a default region size of 100 bp (2 * extendingSize)
#' ## which is extended to include all genomic regions for the closest
#' ## peak to the median position of all peaks included in the region (for each
#' ## experiment).
#' ## A peak from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo(c("chr1", "chr10"), c(249250621, 135534747), NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_FOXA1_01_NarrowPeaks_partial,
#'                         A549_FOSL2_01_NarrowPeaks_partial),
#'     peaks = c(A549_FOXA1_01_Peaks_partial,
#'                         A549_FOSL2_01_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 50,
#'     expandToFitPeakRegion = TRUE,
#'     shrinkToFitPeakRegion = TRUE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL

#' Sites with the greatest evidence of transcription factor binding
#' for the CTCF transcription factor (for demonstration purpose)
#'
#' Sites representing the greatest evidence of enrichment for
#' the CTCF transcription factor (DCC accession: TODO)
#' for regions chr1:246000000-249250621 and chr10:10000000-12500000
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' @name A549_CTCF_MYJ_NarrowPeaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per site.
#'
#' @usage data(A549_CTCF_MYJ_NarrowPeaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_CTCF_MYJ_Peaks_partial}} { the associate
#'                  genomic peaks dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_CTCF_MYJ_NarrowPeaks_partial)
#' data(A549_CTCF_MYJ_Peaks_partial)
#' data(A549_CTCF_MYN_NarrowPeaks_partial)
#' data(A549_CTCF_MYN_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_CTCF_MYJ_Peaks_partial) <- rep("CTCF_MYJ",
#'                               length(A549_CTCF_MYJ_Peaks_partial))
#' names(A549_CTCF_MYJ_NarrowPeaks_partial) <- rep("CTCF_MYJ",
#'                               length(A549_CTCF_MYJ_NarrowPeaks_partial))
#' names(A549_CTCF_MYN_Peaks_partial) <-rep("CTCF_MYN",
#'                               length(A549_CTCF_MYN_Peaks_partial))
#' names(A549_CTCF_MYN_NarrowPeaks_partial) <- rep("CTCF_MYN",
#'                               length(A549_CTCF_MYN_NarrowPeaks_partial))
#'
#' ## Calculating consensus regions for chromosome 10
#' ## with a default region size of 100 bp (2 * extendingSize)
#' ## which is extended to include all genomic regions for the closest
#' ## peak to the median position of all peaks included in the region (for each
#' ## experiment).
#' ## A peak from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo(c("chr10"), c(135534747), NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_CTCF_MYJ_NarrowPeaks_partial,
#'                         A549_CTCF_MYN_NarrowPeaks_partial),
#'     peaks = c(A549_CTCF_MYJ_Peaks_partial,
#'                         A549_CTCF_MYN_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 50,
#'     expandToFitPeakRegion = TRUE,
#'     shrinkToFitPeakRegion = TRUE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL


#' Sites with the greatest evidence of transcription factor binding
#' for the CTCF transcription factor (for demonstration purpose)
#'
#' Sites representing the greatest evidence of enrichment for
#' the CTCF transcription factor (DCC accession: TODO)
#' for regions chr1:246000000-249250621 and chr10:10000000-12500000
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' @name A549_CTCF_MYJ_Peaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per site.
#'
#' @usage data(A549_CTCF_MYJ_Peaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_CTCF_MYJ_NarrowPeaks_partial}} { the associate
#'                  genomic regions dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_CTCF_MYJ_NarrowPeaks_partial)
#' data(A549_CTCF_MYJ_Peaks_partial)
#' data(A549_CTCF_MYN_NarrowPeaks_partial)
#' data(A549_CTCF_MYN_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_CTCF_MYJ_Peaks_partial) <- rep("CTCF_MYJ",
#'                               length(A549_CTCF_MYJ_Peaks_partial))
#' names(A549_CTCF_MYJ_NarrowPeaks_partial) <- rep("CTCF_MYJ",
#'                               length(A549_CTCF_MYJ_NarrowPeaks_partial))
#' names(A549_CTCF_MYN_Peaks_partial) <-rep("CTCF_MYN",
#'                               length(A549_CTCF_MYN_Peaks_partial))
#' names(A549_CTCF_MYN_NarrowPeaks_partial) <- rep("CTCF_MYN",
#'                               length(A549_CTCF_MYN_NarrowPeaks_partial))
#'
#' ## Calculating consensus regions for chromosome 10
#' ## with a default region size of 40 bp (2 * extendingSize)
#' ## which is extended to include all genomic regions for the closest
#' ## peak to the median position of all peaks included in the region (for each
#' ## experiment).
#' ## A peak from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo(c("chr10"), c(135534747), NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_CTCF_MYJ_NarrowPeaks_partial,
#'                         A549_CTCF_MYN_NarrowPeaks_partial),
#'     peaks = c(A549_CTCF_MYJ_Peaks_partial,
#'                         A549_CTCF_MYN_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 20,
#'     expandToFitPeakRegion = FALSE,
#'     shrinkToFitPeakRegion = FALSE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL


#' Sites with the greatest evidence of transcription factor binding
#' for the CTCF transcription factor (for demonstration purpose)
#'
#' Sites representing the greatest evidence of enrichment for
#' the CTCF transcription factor (DCC accession: TODO)
#' for regions chr1:246000000-249250621 and chr10:10000000-12500000
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' @name A549_CTCF_MYN_NarrowPeaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per site.
#'
#' @usage data(A549_CTCF_MYN_NarrowPeaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_CTCF_MYN_Peaks_partial}} { the associate
#'                  genomic regions dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_CTCF_MYJ_NarrowPeaks_partial)
#' data(A549_CTCF_MYJ_Peaks_partial)
#' data(A549_CTCF_MYN_NarrowPeaks_partial)
#' data(A549_CTCF_MYN_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_CTCF_MYJ_Peaks_partial) <- rep("CTCF_MYJ",
#'                               length(A549_CTCF_MYJ_Peaks_partial))
#' names(A549_CTCF_MYJ_NarrowPeaks_partial) <- rep("CTCF_MYJ",
#'                               length(A549_CTCF_MYJ_NarrowPeaks_partial))
#' names(A549_CTCF_MYN_Peaks_partial) <-rep("CTCF_MYN",
#'                               length(A549_CTCF_MYN_Peaks_partial))
#' names(A549_CTCF_MYN_NarrowPeaks_partial) <- rep("CTCF_MYN",
#'                               length(A549_CTCF_MYN_NarrowPeaks_partial))
#'
#' ## Calculating consensus regions for chromosome 10
#' ## with a default region size of 40 bp (2 * extendingSize)
#' ## which is extended to include all genomic regions for the closest
#' ## peak to the median position of all peaks included in the region (for each
#' ## experiment).
#' ## A peak from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo("chr10", 135534747, NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_CTCF_MYJ_NarrowPeaks_partial,
#'                         A549_CTCF_MYN_NarrowPeaks_partial),
#'     peaks = c(A549_CTCF_MYJ_Peaks_partial,
#'                         A549_CTCF_MYN_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 20,
#'     expandToFitPeakRegion = FALSE,
#'     shrinkToFitPeakRegion = FALSE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL


#' Sites with the greatest evidence of transcription factor binding
#' for the CTCF transcription factor (for demonstration purpose)
#'
#' Sites representing the greatest evidence of enrichment for
#' the CTCF transcription factor (DCC accession: TODO)
#' for regions chr1:246000000-249250621 and chr10:10000000-12500000
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' @name A549_CTCF_MYN_Peaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per site.
#'
#' @usage data(A549_CTCF_MYN_Peaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_CTCF_MYN_NarrowPeaks_partial}} { the associate
#'                  genomic regions dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_CTCF_MYJ_NarrowPeaks_partial)
#' data(A549_CTCF_MYJ_Peaks_partial)
#' data(A549_CTCF_MYN_NarrowPeaks_partial)
#' data(A549_CTCF_MYN_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_CTCF_MYJ_Peaks_partial) <- rep("CTCF_MYJ",
#'                               length(A549_CTCF_MYJ_Peaks_partial))
#' names(A549_CTCF_MYJ_NarrowPeaks_partial) <- rep("CTCF_MYJ",
#'                               length(A549_CTCF_MYJ_NarrowPeaks_partial))
#' names(A549_CTCF_MYN_Peaks_partial) <-rep("CTCF_MYN",
#'                               length(A549_CTCF_MYN_Peaks_partial))
#' names(A549_CTCF_MYN_NarrowPeaks_partial) <- rep("CTCF_MYN",
#'                               length(A549_CTCF_MYN_NarrowPeaks_partial))
#'
#' ## Calculating consensus regions for chromosomes 1
#' ## with a default region size of 40 bp (2 * extendingSize)
#' ## which is extended to include all genomic regions for the closest
#' ## peak to the median position of all peaks included in the region (for each
#' ## experiment).
#' ## A peak from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo(c("chr1"), c(249250621), NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_CTCF_MYJ_NarrowPeaks_partial,
#'                         A549_CTCF_MYN_NarrowPeaks_partial),
#'     peaks = c(A549_CTCF_MYJ_Peaks_partial,
#'                         A549_CTCF_MYN_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 20,
#'     expandToFitPeakRegion = FALSE,
#'     shrinkToFitPeakRegion = FALSE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL


#' Sites with the greatest evidence of transcription factor binding
#' for the NR3C1 transcription factor from ENCODE (DDC accession: ENCFF002CFQ).
#' For demonstration purpose.
#'
#' Sites representing the greatest evidence of enrichment for
#' the NR3C1 transcription factor (DCC accession: ENCFF002CFQ)
#' for regions chr2:40000000-50000000 and chr3:10000000-13000000
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' The peaks and ranges have been obtained using an optimal IDR analysis
#' done on all replicates.
#'
#' @name A549_NR3C1_CFQ_Peaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per site. The peaks are
#' surronded by ranges present in the dataset
#' \code{A549_NR3C1_CFQ_NarrowPeaks_partial}.
#'
#' @usage data(A549_NR3C1_CFQ_Peaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_NR3C1_CFQ_NarrowPeaks_partial}} { the associate
#'                  genomic regions dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_NR3C1_CFQ_NarrowPeaks_partial)
#' data(A549_NR3C1_CFQ_Peaks_partial)
#' data(A549_NR3C1_CFS_NarrowPeaks_partial)
#' data(A549_NR3C1_CFS_Peaks_partial)
#' data(A549_NR3C1_CFR_NarrowPeaks_partial)
#' data(A549_NR3C1_CFR_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_NR3C1_CFQ_NarrowPeaks_partial) <- rep("NR3C1_CFQ",
#'                               length(A549_NR3C1_CFQ_NarrowPeaks_partial))
#' names(A549_NR3C1_CFQ_Peaks_partial) <- rep("NR3C1_CFQ",
#'                               length(A549_NR3C1_CFQ_Peaks_partial))
#' names(A549_NR3C1_CFS_NarrowPeaks_partial) <-rep("NR3C1_CFS",
#'                               length(A549_NR3C1_CFS_NarrowPeaks_partial))
#' names(A549_NR3C1_CFS_Peaks_partial) <- rep("NR3C1_CFS",
#'                               length(A549_NR3C1_CFS_Peaks_partial))
#' names(A549_NR3C1_CFR_NarrowPeaks_partial) <-rep("NR3C1_CFR",
#'                               length(A549_NR3C1_CFR_NarrowPeaks_partial))
#' names(A549_NR3C1_CFR_Peaks_partial) <- rep("NR3C1_CFR",
#'                               length(A549_NR3C1_CFR_Peaks_partial))
#'
#' ## Calculating consensus regions for chromosome 3
#' ## with a default region size of 140 bp (2 * extendingSize)
#' ## which is extended to include all genomic regions for the closest
#' ## peak to the median position of all peaks included in the region (for
#' ## each experiment).
#' ## Peaks from at least 2 experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo(c("chr3"), c(198022430), NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_NR3C1_CFQ_NarrowPeaks_partial,
#'                         A549_NR3C1_CFS_NarrowPeaks_partial,
#'                         A549_NR3C1_CFR_NarrowPeaks_partial),
#'     peaks = c(A549_NR3C1_CFQ_Peaks_partial,
#'                         A549_NR3C1_CFS_Peaks_partial,
#'                         A549_NR3C1_CFR_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 70,
#'     expandToFitPeakRegion = FALSE,
#'     shrinkToFitPeakRegion = FALSE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL


#' Ranges with the greatest evidence of transcription factor binding
#' for the NR3C1 transcription factor from ENCODE (DDC accession: ENCFF002CFQ).
#' For demonstration purpose.
#'
#' Ranges representing the greatest evidence of enrichment for
#' the NR3C1 transcription factor (DCC accession: ENCFF002CFQ)
#' for regions chr2:40000000-50000000 and chr3:10000000-13000000
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' The peaks and ranges have been obtained using an optimal IDR analysis
#' done on all replicates.
#'
#' @name A549_NR3C1_CFQ_NarrowPeaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per site. The ranges are
#' surronding the peaks present in the dataset
#' \code{A549_NR3C1_CFQ_Peaks_partial}.
#'
#' @usage data(A549_NR3C1_CFQ_NarrowPeaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_NR3C1_CFQ_Peaks_partial}} { the associate
#'                  genomic regions dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_NR3C1_CFQ_NarrowPeaks_partial)
#' data(A549_NR3C1_CFQ_Peaks_partial)
#' data(A549_NR3C1_CFS_NarrowPeaks_partial)
#' data(A549_NR3C1_CFS_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_NR3C1_CFQ_NarrowPeaks_partial) <- rep("NR3C1_CFQ",
#'                               length(A549_NR3C1_CFQ_NarrowPeaks_partial))
#' names(A549_NR3C1_CFQ_Peaks_partial) <- rep("NR3C1_CFQ",
#'                               length(A549_NR3C1_CFQ_Peaks_partial))
#' names(A549_NR3C1_CFS_NarrowPeaks_partial) <-rep("NR3C1_CFS",
#'                               length(A549_NR3C1_CFS_NarrowPeaks_partial))
#' names(A549_NR3C1_CFS_Peaks_partial) <- rep("NR3C1_CFS",
#'                               length(A549_NR3C1_CFS_Peaks_partial))
#'
#' ## Calculating consensus regions for chromosome 3
#' ## with a default region size of 300 bp (2 * extendingSize)
#' ## which is extended to include all genomic regions for the closest
#' ## peak to the median position of all peaks included in the region (for
#' ## each experiment).
#' ## Peaks from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo(c("chr3"), c(198022430), NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_NR3C1_CFQ_NarrowPeaks_partial,
#'                         A549_NR3C1_CFS_NarrowPeaks_partial),
#'     peaks = c(A549_NR3C1_CFQ_Peaks_partial,
#'                         A549_NR3C1_CFS_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 150,
#'     expandToFitPeakRegion = FALSE,
#'     shrinkToFitPeakRegion = TRUE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL


#' Sites with the greatest evidence of transcription factor binding
#' for the NR3C1 transcription factor from ENCODE (DDC accession: ENCFF002CFR).
#' For demonstration purpose.
#'
#' Sites representing the greatest evidence of enrichment for
#' the NR3C1 transcription factor (DCC accession: ENCFF002CFR)
#' for regions chr2:40000000-50000000 and chr3:10000000-13000000
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' The peaks and ranges have been obtained using an optimal IDR analysis
#' done on all replicates.
#'
#' @name A549_NR3C1_CFR_Peaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per site. The peaks are
#' surronded by ranges present in the dataset
#' \code{A549_NR3C1_CFR_NarrowPeaks_partial}.
#'
#' @usage data(A549_NR3C1_CFR_Peaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_NR3C1_CFR_NarrowPeaks_partial}} { the associate
#'                  genomic regions dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_NR3C1_CFQ_NarrowPeaks_partial)
#' data(A549_NR3C1_CFQ_Peaks_partial)
#' data(A549_NR3C1_CFR_NarrowPeaks_partial)
#' data(A549_NR3C1_CFR_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_NR3C1_CFQ_NarrowPeaks_partial) <- rep("NR3C1_CFQ",
#'                               length(A549_NR3C1_CFQ_NarrowPeaks_partial))
#' names(A549_NR3C1_CFQ_Peaks_partial) <- rep("NR3C1_CFQ",
#'                               length(A549_NR3C1_CFQ_Peaks_partial))
#' names(A549_NR3C1_CFR_NarrowPeaks_partial) <-rep("NR3C1_CFR",
#'                               length(A549_NR3C1_CFR_NarrowPeaks_partial))
#' names(A549_NR3C1_CFR_Peaks_partial) <- rep("NR3C1_CFR",
#'                               length(A549_NR3C1_CFR_Peaks_partial))
#'
#' ## Calculating consensus regions for chromosome 2
#' ## with a default region size of 40 bp (2 * extendingSize)
#' ## which is extended to include all genomic regions for the closest
#' ## peak to the median position of all peaks included in the region (for
#' ## each experiment).
#' ## Peaks from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo(c("chr2"), c(243199373), NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_NR3C1_CFQ_NarrowPeaks_partial,
#'                         A549_NR3C1_CFR_NarrowPeaks_partial),
#'     peaks = c(A549_NR3C1_CFQ_Peaks_partial,
#'                         A549_NR3C1_CFR_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 20,
#'     expandToFitPeakRegion = TRUE,
#'     shrinkToFitPeakRegion = FALSE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL


#' Ranges with the greatest evidence of transcription factor binding
#' for the NR3C1 transcription factor from ENCODE (DDC accession: ENCFF002CFR).
#' For demonstration purpose.
#'
#' Ranges representing the greatest evidence of enrichment for
#' the NR3C1 transcription factor (DCC accession: ENCFF002CFR)
#' for regions chr2:40000000-50000000 and chr3:10000000-13000000
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' The peaks and ranges have been obtained using an optimal IDR analysis
#' done on all replicates.
#'
#' @name A549_NR3C1_CFR_NarrowPeaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per site. The ranges are
#' surronding the peaks present in the dataset
#' \code{A549_NR3C1_CFR_Peaks_partial}.
#'
#' @usage data(A549_NR3C1_CFR_NarrowPeaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_NR3C1_CFR_Peaks_partial}} { the associate
#'                  genomic regions dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_NR3C1_CFQ_NarrowPeaks_partial)
#' data(A549_NR3C1_CFQ_Peaks_partial)
#' data(A549_NR3C1_CFR_NarrowPeaks_partial)
#' data(A549_NR3C1_CFR_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_NR3C1_CFQ_NarrowPeaks_partial) <- rep("NR3C1_CFQ",
#'                               length(A549_NR3C1_CFQ_NarrowPeaks_partial))
#' names(A549_NR3C1_CFQ_Peaks_partial) <- rep("NR3C1_CFQ",
#'                               length(A549_NR3C1_CFQ_Peaks_partial))
#' names(A549_NR3C1_CFR_NarrowPeaks_partial) <-rep("NR3C1_CFR",
#'                               length(A549_NR3C1_CFR_NarrowPeaks_partial))
#' names(A549_NR3C1_CFR_Peaks_partial) <- rep("NR3C1_CFR",
#'                               length(A549_NR3C1_CFR_Peaks_partial))
#'
#' ## Calculating consensus regions for chromosome 2
#' ## with a default region size of 250 bp (2 * extendingSize)
#' ## which is extended to include all genomic regions for the closest
#' ## peak to the median position of all peaks included in the region (for
#' ## each experiment).
#' ## Peaks from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo(c("chr2"), c(243199373), NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_NR3C1_CFQ_NarrowPeaks_partial,
#'                         A549_NR3C1_CFR_NarrowPeaks_partial),
#'     peaks = c(A549_NR3C1_CFQ_Peaks_partial,
#'                         A549_NR3C1_CFR_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 125,
#'     expandToFitPeakRegion = TRUE,
#'     shrinkToFitPeakRegion = TRUE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL

#' Sites with the greatest evidence of transcription factor binding
#' for the NR3C1 transcription factor from ENCODE (DDC accession: ENCFF002CFS).
#' For demonstration purpose.
#'
#' Sites representing the greatest evidence of enrichment for
#' the NR3C1 transcription factor (DCC accession: ENCFF002CFS)
#' for regions chr2:40000000-50000000 and chr3:10000000-13000000
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' The peaks and ranges have been obtained using an optimal IDR analysis
#' done on all replicates.
#'
#' @name A549_NR3C1_CFS_Peaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per site. The peaks are
#' surronded by ranges present in the dataset
#' \code{A549_NR3C1_CFS_NarrowPeaks_partial}.
#'
#' @usage data(A549_NR3C1_CFS_Peaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_NR3C1_CFS_NarrowPeaks_partial}} { the associate
#'                  genomic regions dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_NR3C1_CFQ_NarrowPeaks_partial)
#' data(A549_NR3C1_CFQ_Peaks_partial)
#' data(A549_NR3C1_CFS_NarrowPeaks_partial)
#' data(A549_NR3C1_CFS_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_NR3C1_CFQ_NarrowPeaks_partial) <- rep("NR3C1_CFQ",
#'                               length(A549_NR3C1_CFQ_NarrowPeaks_partial))
#' names(A549_NR3C1_CFQ_Peaks_partial) <- rep("NR3C1_CFQ",
#'                               length(A549_NR3C1_CFQ_Peaks_partial))
#' names(A549_NR3C1_CFS_NarrowPeaks_partial) <-rep("NR3C1_CFS",
#'                               length(A549_NR3C1_CFS_NarrowPeaks_partial))
#' names(A549_NR3C1_CFS_Peaks_partial) <- rep("NR3C1_CFS",
#'                               length(A549_NR3C1_CFS_Peaks_partial))
#'
#' ## Calculating consensus regions for chromosome 2
#' ## with a default region size of 80 bp (2 * extendingSize)
#' ## which is extended to include all genomic regions for the closest
#' ## peak to the median position of all peaks included in the region (for
#' ## each experiment).
#' ## Peaks from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo(c("chr2"), c(243199373), NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_NR3C1_CFQ_NarrowPeaks_partial,
#'                         A549_NR3C1_CFS_NarrowPeaks_partial),
#'     peaks = c(A549_NR3C1_CFQ_Peaks_partial,
#'                         A549_NR3C1_CFS_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 40,
#'     expandToFitPeakRegion = FALSE,
#'     shrinkToFitPeakRegion = FALSE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL


#' Ranges with the greatest evidence of transcription factor binding
#' for the NR3C1 transcription factor from ENCODE (DDC accession: ENCFF002CFS).
#' For demonstration purpose.
#'
#' Ranges representing the greatest evidence of enrichment for
#' the NR3C1 transcription factor (DCC accession: ENCFF002CFS)
#' for regions chr2:40000000-50000000 and chr3:10000000-13000000
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' The peaks and ranges have been obtained using an optimal IDR analysis
#' done on all replicates.
#'
#' @name A549_NR3C1_CFS_NarrowPeaks_partial
#'
#' @docType data
#'
#' @format A \code{GRanges} containing one entry per site. The ranges are
#' surronding the peaks present in the dataset
#' \code{A549_NR3C1_CFs_Peaks_partial}.
#'
#' @usage data(A549_NR3C1_CFS_NarrowPeaks_partial)
#'
#' @references
#' \itemize{
#'     \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia
#' of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{A549_NR3C1_CFS_Peaks_partial}} { the associate
#'                  genomic regions dataset.}
#'     \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading datasets
#' data(A549_NR3C1_CFQ_NarrowPeaks_partial)
#' data(A549_NR3C1_CFQ_Peaks_partial)
#' data(A549_NR3C1_CFS_NarrowPeaks_partial)
#' data(A549_NR3C1_CFS_Peaks_partial)
#'
#' ## Assigning experiment name to each row of the dataset.
#' ## NarrowPeak and Peak datasets from the same experiment must
#' ## have identical names.
#' names(A549_NR3C1_CFQ_NarrowPeaks_partial) <- rep("NR3C1_CFQ",
#'                               length(A549_NR3C1_CFQ_NarrowPeaks_partial))
#' names(A549_NR3C1_CFQ_Peaks_partial) <- rep("NR3C1_CFQ",
#'                               length(A549_NR3C1_CFQ_Peaks_partial))
#' names(A549_NR3C1_CFS_NarrowPeaks_partial) <-rep("NR3C1_CFS",
#'                               length(A549_NR3C1_CFS_NarrowPeaks_partial))
#' names(A549_NR3C1_CFS_Peaks_partial) <- rep("NR3C1_CFS",
#'                               length(A549_NR3C1_CFS_Peaks_partial))
#'
#' ## Calculating consensus regions for chromosome 2
#' ## with a default region size of 500 bp (2 * extendingSize)
#' ## which is extended to include all genomic regions for the closest
#' ## peak to the median position of all peaks included in the region (for
#' ## each experiment).
#' ## Peaks from both experiments must be present in a region to
#' ## be retained as a consensus region.
#' chrList <- Seqinfo(c("chr2"), c(243199373), NA)
#' findConsensusPeakRegions(
#'     narrowPeaks = c(A549_NR3C1_CFQ_NarrowPeaks_partial,
#'                         A549_NR3C1_CFS_NarrowPeaks_partial),
#'     peaks = c(A549_NR3C1_CFQ_Peaks_partial,
#'                         A549_NR3C1_CFS_Peaks_partial),
#'     chrInfo = chrList,
#'     extendingSize = 500,
#'     expandToFitPeakRegion = FALSE,
#'     shrinkToFitPeakRegion = FALSE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
NULL