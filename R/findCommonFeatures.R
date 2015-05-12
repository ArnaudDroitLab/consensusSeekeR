#' @title Extract regions sharing the same features in more than one experiment
#' 
#' @description Find regions sharing the same features for a minimum number of
#'          experiments using narrowPeak files. The peaks and narrow regions
#'          are extracted from the narrowPeaks files and used to identify 
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
#' @param narrowPeakFiles a \code{vector} containing the narrowPeak files to
#'          use.
#' @param chrList a \code{vector} containing the name of the chromosomes to 
#'          analyze or the name \code{"ALL"} which indicate that all
#'          chromosomes must be analyzed. When \code{NULL}, no
#'          new term is added. Default : \code{NULL}.
#' @param extendingSize a \code{numeric} value indicating the size of padding 
#'          at each side of the peaks median position to create the consensus
#'          region. The minimum size of the consensu region will be equal to
#'          twice the value of the \code{extendingSize} parameter. The size of 
#'          the \code{extendingSize} must be a positive integer. Default = 250.
#' @param includeAllPeakRegion a \code{logical} indicating if the region size,
#'          which is set by the \code{extendingSize} parameter is extended to 
#'          include all region of the peak closest to the peaks median 
#'          position for each experiment.
#' @param minNbrExp a \code{numeric} indicating the minimum number of BED files
#'          in which a peak must be present for a region to be retained. The
#'          numeric must be a positive integer inferior or equal to the number 
#'          of files present in the \code{narrowpeaksBEDFiles} parameter.
#'          Default = 1.
#' @param nbrThreads a \code{numeric} indicating the number of threads to use
#'          in parallel. The \code{nbrThreads} must be a positive integer. 
#'          Default = 1.
#' 
#' @return an object of \code{class} "commonFeatures". 
#' 
#' @author Astrid Louise Deschenes
#' @importFrom BiocGenerics start end
#' @importFrom stringr str_split
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam 
#'                  multicoreWorkers
#' @export
findConsensusPeakRegions <- function(narrowPeakFiles, chrList = "ALL", 
                               extendingSize = 250, 
                               includeAllPeakRegion = TRUE, minNbrExp = 1, 
                               nbrThreads = 1) {
    
    # Parameters validation
    findConsensusPeakRegionsValidation(narrowPeakFiles, chrList, extendingSize,
                                includeAllPeakRegion, minNbrExp, nbrThreads)
    
    # Create objects that are going to contain the final extracted values
    allPeaks <- GRanges()
    allNarrowPeaks <- GRanges()
    
    # Extract peaks and regions from each file present in the file vector
    for (files in narrowPeakFiles) {
        data <- readNarrowPeak(files)
        allPeaks <- append(allPeaks, data$peak)
        allNarrowPeaks <- append(allNarrowPeaks, data$narrowPeak)
    }
    
    # Select the type of object used for parallel processing
    coreParam <- MulticoreParam(workers = nbrThreads)
    if (nbrThreads == 1 || multicoreWorkers() == 1) {
        coreParam <- SerialParam()
    }
    
    # Extract the list of chromosomes to analyse
    allChr <- levels(seqnames(allPeaks))
    if (chrList == "ALL") {
        # The list of chromosomes correspond to the global list
        chrList <- allChr 
    } else {
        # The list of chromosomes correspond to the chromosomes from the
        # specified parameter which are present in the data
        chrList <- subset(allChr, allChr %in% chrList)
    }
    
    # At least one chromosome must be analyzed
    if (length(chrList) == 0) {
        stop("No chromosome correspond to the given list: ", chrList)
    }
    
    # Process to regions extraction using parallel threads when available
    results <- bplapply(chrList, 
                FUN = findConsensusPeakRegionsForOneChrom,
                allPeaks = allPeaks, allNarrowPeaks = allNarrowPeaks, 
                extendingSize = extendingSize, includeAllPeakRegion =
                includeAllPeakRegion, minNbrExp = minNbrExp, 
                BPPARAM = coreParam)
    
    # Merge extracted regions
    finalRegions <- GRanges()
    for (i in 1:length(results)) {
        finalRegions<-c(finalRegions, results[[i]][["features"]])
    }
    
    return(finalRegions)
}