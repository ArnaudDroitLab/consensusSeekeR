#' @title Extract regions sharing the same features in more than one experiment
#'
#' @description Find regions sharing the same features for a minimum number of
#'          experiments using called peaks of signal enrichment based on
#'          pooled, normalized data (mainly coming from narrowPeak files).
#'          The peaks and narrow peaks are used to identify
#'          the consensus regions. The minimum number of experiments that must
#'          have at least on peak in a region so that it is retained as a
#'          consensus region is specified by user, as well as the size of
#'          mining regions. Only the
#'          chromosomes specified by the user are treated. The function can be
#'          parallized by specifying a number of threads superior to 1.
#'          However, Windows does not support multicore  evaluation.
#'
#'          When the padding is small, the detected regions are smaller than
#'          the one that could be obtained by doing an overlap of the narrow
#'          regions. Even more, the parameter specifying the minimum number of
#'          experiments needed to retain a region add versatility to the
#'          function.
#'
#'          Beware that the side of the padding can have a large effect on
#'          the detected consensus
#'          regions. It is recommanded to test more than one size and to do
#'          some manual validation of the resulting consensus regions before
#'          selecting the final padding size.
#'
#' @param narrowPeaks a \code{GRanges} containing
#'          called peak regions of signal enrichment based on pooled,
#'          normalized data
#'          for all analyzed experiments. All \code{GRanges} entries must
#'          have a metadata field called "name" which identifies the region to
#'          the called peak. All \code{GRanges} entries must also
#'          have a row name which identifies the experiment of origin. Each
#'          \code{peaks} entry must have
#'          an associated \code{narrowPeaks} entry. A \code{GRanges} entry is
#'          associated to a \code{narrowPeaks} entry by having a identical
#'          metadata "name" field and a identical row name.
#' @param peaks a \code{GRanges} containing
#'          called peaks of signal enrichment based on pooled,
#'          normalized data
#'          for all analyzed experiments. All \code{GRanges} entries must
#'          have a metadata field called "name" which identifies the called
#'          peak. All \code{GRanges} entries must
#'          have a row name which identifies the experiment of origin. Each
#'          \code{peaks} entry must have
#'          an associated \code{narrowPeaks} entry. A \code{GRanges} entry is
#'          associated to a \code{narrowPeaks} entry by having a identical
#'          metadata "name" field and a identical row name.
#' @param chrInfo a \code{Seqinfo} containing the name and the length of the
#'          chromosomes to analyze. Only the chomosomes contained in this
#'          \code{Seqinfo} will be analyzed.
#' @param extendingSize a \code{numeric} value indicating the size of padding
#'          on both sides of the position of the peaks median to create the
#'          consensus region. The minimum size of the consensus region is
#'          equal to
#'          twice the value of the \code{extendingSize} parameter. The size of
#'          the \code{extendingSize} must be a positive integer. Default = 250.
#' @param includeAllPeakRegion a \code{logical} indicating if the region size,
#'          which is set by the \code{extendingSize} parameter is extended to
#'          include the entire narrow peak region of the peak closest to the
#'          position of the peaks median
#'          for each experiment. When two peaks are at equidistance
#'          of the peaks median for a specific experiment, both peaks are
#'          used to extend the final consensus region.
#' @param minNbrExp a \code{numeric} or a \code{integer} indicating the
#'          minimum number of experiments
#'          in which at least one peak must be present for a potential
#'          consensus region. The
#'          numeric must be a positive integer inferior or equal to the number
#'          of experiments present in the \code{narrowPeaks} and \code{peaks}
#'          parameters.
#'          Default = 1.
#' @param nbrThreads a \code{numeric} or a \code{integer} indicating the
#'          number of threads to use
#'          in parallel. The \code{nbrThreads} must be a positive integer.
#'          Default = 1.
#'
#' @return an object of \code{class} "consensusRanges" containing :
#'      \itemize{
#'          \item \code{call} the matched call
#'          \item \code{consensusRanges} a \code{GRanges} containing the
#'                      consensus regions.
#'      }
#'
#' @author Astrid Louise Deschenes
#' @importFrom BiocGenerics start end
#' @importFrom stringr str_split
#' @importFrom IRanges IRanges unlist
#' @importFrom GenomicRanges GRanges GRangesList .__T__split:base
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam
#'                  multicoreWorkers bpmapply
#' @importFrom GenomeInfoDb Seqinfo seqinfo seqnames
#' @export
findConsensusPeakRegions <- function(narrowPeaks, peaks, chrInfo,
                               extendingSize = 250,
                               includeAllPeakRegion = TRUE, minNbrExp = 1,
                               nbrThreads = 1) {
    # Get call information
    cl <- match.call()

    # Parameters validation
    findConsensusPeakRegionsValidation(narrowPeaks, peaks, chrInfo,
            extendingSize, includeAllPeakRegion, minNbrExp, nbrThreads)

    # Select the type of object used for parallel processing
    coreParam <- MulticoreParam(workers = nbrThreads)
    if (nbrThreads == 1 || multicoreWorkers() == 1) {
        coreParam <- SerialParam()
    }

    # Preparing data
    narrowPeaksSplit <- GenomicRanges::split(narrowPeaks,
                                                seqnames(narrowPeaks))
    peaksSplit <- GenomicRanges::split(peaks, seqnames(peaks))
    rm(peaks)
    rm(narrowPeaks)
    selectedNarrowPeaksSplit <- narrowPeaksSplit[names(narrowPeaksSplit) %in%
                                                seqnames(chrInfo)]
    selectedPeaksSplit <- peaksSplit[names(peaksSplit) %in%
                                                seqnames(chrInfo)]
    rm(narrowPeaksSplit)
    rm(peaksSplit)

    # Running each chromosome on a separate thread
    results2 <- bpmapply(findConsensusPeakRegionsForOneChrom,
                        chrName = seqnames(chrInfo),
                        allPeaks = selectedPeaksSplit,
                        allNarrowPeaks = selectedNarrowPeaksSplit,
                        MoreArgs = c(extendingSize = extendingSize,
                        includeAllPeakRegion =
                        includeAllPeakRegion, minNbrExp = minNbrExp,
                        chrList = chrInfo),
                        BPPARAM = coreParam)

    # Creating result list
    z <- list(call = cl,
                    consensusRanges = IRanges::unlist(GRangesList((results2)),
                    recursive = TRUE, use.names = FALSE))

    class(z)<-"consensusRanges"

    return(z)
}
