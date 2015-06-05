#' @title Extract regions sharing the same features in more than one experiment
#' 
#' @description Find regions sharing the same features for a minimum number of
#'          experiments using called peaks of signal enrichment based on 
#'          pooled, normalized data (mainly coming from narrowPeak files). 
#'          The peaks and narrow peaks and used to identify 
#'          the common regions. The minimum number of experiments that must a
#'          peak in a common regions for that region to be selected is 
#'          specified by user, as well as the size of padding. Only the 
#'          chromosomes specified by the user are treated. The function can be 
#'          parallized by specifying a number of threads superior to 1. 
#'          
#'          When the padding is small, the detected regions are smaller than 
#'          the one that could be obtained by doing an overlap of the narrow
#'          regions. Even more, the parameter specifying the minimum number of 
#'          experiments needed to retain a region add versatility to the 
#'          function.
#'          
#'          The side of the padding can have a large effect on the detected
#'          regions. It is recommanded to test more than one size and to do
#'          some manual validation of the resulting regions before selecting
#'          the final padding size.
#' 
#' @param narrowPeaks a \code{vector} containing \code{GRanges} representing 
#'          called peaks of signal enrichment based on pooled, normalized data 
#'          for all experiments.
#' @param peaks a \code{vector} containing \code{GRanges} representing peaks.
#' @param chrInfo a \code{Seqinfo} containing the name and the length of the 
#'          chromosomes to analyze.
#' @param extendingSize a \code{numeric} value indicating the size of padding 
#'          at each side of the peaks median position to create the consensus
#'          region. The minimum size of the consensu region will be equal to
#'          twice the value of the \code{extendingSize} parameter. The size of 
#'          the \code{extendingSize} must be a positive integer. Default = 250.
#' @param includeAllPeakRegion a \code{logical} indicating if the region size,
#'          which is set by the \code{extendingSize} parameter is extended to 
#'          include all region of the peak closest to the peaks median 
#'          position for each experiment. When two peaks are at equal distance
#'          of the peaks median for one experiment, both peaks are used to
#'          extend the final consensus region.
#' @param minNbrExp a \code{numeric} indicating the minimum number of BED files
#'          in which a peak must be present for a region to be retained. The
#'          numeric must be a positive integer inferior or equal to the number 
#'          of files present in the \code{narrowPeakFiles} parameter.
#'          Default = 1.
#' @param nbrThreads a \code{numeric} indicating the number of threads to use
#'          in parallel. The \code{nbrThreads} must be a positive integer. 
#'          Default = 1.
#' 
#' @return an object of \code{class} "consensusRanges". 
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
    cl <- match.call()
    
    # Parameters validation
    findConsensusPeakRegionsValidation(narrowPeaks, peaks, chrInfo, 
            extendingSize, includeAllPeakRegion, minNbrExp, nbrThreads)
    
    # Select the type of object used for parallel processing
    coreParam <- MulticoreParam(workers = nbrThreads)
    if (nbrThreads == 1 || multicoreWorkers() == 1) {
        coreParam <- SerialParam()
    }
    
    # Process to regions extraction using parallel threads when available
#     results <- bplapply(names(chrInfo), 
#                 FUN = findConsensusPeakRegionsForOneChrom,
#                 allPeaks = peaks, allNarrowPeaks = narrowPeaks, 
#                 extendingSize = extendingSize, includeAllPeakRegion =
#                 includeAllPeakRegion, minNbrExp = minNbrExp, chrList = chrInfo,
#                 BPPARAM = coreParam)
    
    narrowPeaksSplit <- GenomicRanges::split(narrowPeaks, seqnames(narrowPeaks))
    peaksSplit <- GenomicRanges::split(peaks, seqnames(peaks))
    rm(peaks)
    rm(narrowPeaks)
    selectedNarrowPeaksSplit <- narrowPeaksSplit[names(narrowPeaksSplit) %in% 
                                                    seqnames(chrInfo)]
    selectedPeaksSplit <- peaksSplit[names(peaksSplit) %in% 
                                                 seqnames(chrInfo)]
    rm(narrowPeaksSplit)
    rm(peaksSplit)
    
    results2 <- bpmapply(findConsensusPeakRegionsForOneChrom,
                        chrName = seqnames(chrInfo),
                        allPeaks = selectedPeaksSplit, 
                        allNarrowPeaks = selectedNarrowPeaksSplit, 
                        MoreArgs = c(extendingSize = extendingSize, 
                        includeAllPeakRegion =
                        includeAllPeakRegion, minNbrExp = minNbrExp,
                        chrList = chrInfo),
                        BPPARAM = coreParam)
    
                             
    z <- list(call = cl,
                    consensusRanges = IRanges::unlist(GRangesList((results2)), 
                    recursive = TRUE, use.names = FALSE))


    class(z)<-"consensusRanges"

    return(z)
}
