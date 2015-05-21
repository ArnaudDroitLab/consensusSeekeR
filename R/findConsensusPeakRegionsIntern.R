#' @title Parameter validation for the \code{\link{findConsensusPeakRegions}} 
#'      function
#' 
#' @description Validation of all parameters needed by the public
#'      \code{\link{findConsensusPeakRegions}} function.
#' 
#' @param narrowPeakFiles a \code{vector} containing the narrowPeak files to
#'          use for the regions selection.
#' @param chrList a \code{vector} containing the name of the chromosomes to 
#'          analyze or the name \code{"ALL"} which indicate that all
#'          chromosomes must be analyzed. When \code{NULL}, no
#'          new term is added. Default : \code{NULL}.
#' @param extendingSize a \code{numeric} value indicating the size of padding 
#'          at each side of the peaks median position to create the consensus
#'          region. The minimum size of the consensu region will be equal to
#'          twice the value of the \code{extendingSize} parameter. The size of 
#'          the \code{extendingSize} must be a positive integer. Default = 250.
#' @param includeAllPeakRegion a \code{logical} indicating if the region set by
#'          the \code{extendingSize} parameter is extended to include all
#'          region of the peak closest to the peaks median position for each
#'          experiment.
#' @param minNbrExp a \code{numeric} indicating the minimum number of BED files
#'          in which a peak must be present for a region to be retained. The
#'          numeric must be a positive value inferior or equal to the number of 
#'          files present in the \code{narrowPeakFiles} parameter.
#'          Default = 1.
#' @param nbrThreads a \code{numeric} indicating the number of threads to use
#'          in parallel.
#' 
#' @return \code{0} indicating that all parameters validations have been
#'      successful.
#' 
#' @author Astrid Louise Deschenes
#' @keywords internal
findConsensusPeakRegionsValidation <- function(narrowPeakFiles, chrList, 
                                        extendingSize, includeAllPeakRegion, 
                                        minNbrExp, nbrThreads) {
    
    if (is.vector(narrowPeakFiles) && (!is.character(narrowPeakFiles) 
            || !all(sapply(narrowPeakFiles, file.exists)))) {
        stop("peaksBEDlist must be a vector of existing BED files")
    }
    
    if (chrList != "ALL" && !(is.vector(chrList) && is.character(chrList))) {
        stop(paste0("chrList must either be the value \"ALL\" or a ",
             "vector of chromosomes names"))
    }
    
    if (!isInteger(extendingSize) || extendingSize < 1 ) {
        stop("extendingSize must be a non-negative integer")
    }
    
    if (!is.logical(includeAllPeakRegion)) {
        stop("includeAllPeakRegion must be a logical value")
    }
    
    if (!isInteger(minNbrExp) || minNbrExp < 1  || 
            minNbrExp > length(narrowPeakFiles)) {
        stop(paste0("minNbrExp must be a non-negative integer inferior or ", 
                "equal to the number of experiments."))
    }
    
    if (!isInteger(nbrThreads) || nbrThreads < 1 ) {
        stop("nbrThreads must be a non-negative integer")
    }
    
    return(0)
}

#' @title Validate if a value is an integer
#' 
#' @description Validate if the value passed to the function is an integer or 
#'          not.
#'
#' @param value an object to validate.
#' 
#' @return \code{TRUE} is the parameter is a integer; otherwise \code{FALSE}
#' 
#' @author Astrid Louise Deschenes
#' @keywords internal
isInteger <- function(value) {
    return(is.integer(value) || (is.numeric(value) && 
                as.integer(value) == value))
}

#' @title TODO
#' 
#' @description TODO
#' 
#' @param chrName the name of the chromosome to analyse.
#' @param extendingSize a \code{numeric} value indicating the size of padding 
#'          at each side of the peaks median position to create the consensus
#'          region. The minimum size of the consensu region will be equal to
#'          twice the value of the \code{extendingSize} parameter. The size of 
#'          the \code{extendingSize} must be a positive integer. Default = 250.
#' @param includeAllPeakRegion a \code{logical} indicating if the region set by
#'          the \code{extendingSize} parameter is extended to include all
#'          region of the peak closest to the peaks median position for each
#'          experiment.
#' @param minNbrExp a \code{numeric} indicating the minimum number of BED files
#'          in which a peak must be present for a region to be retained. The
#'          numeric must be a positive value inferior or equal to the number of 
#'          files present in the \code{narrowPeakFiles} parameter.
#'          Default = 1.
#' @param allPeaks a \code{GRanges} containing all peaks from all experiments
#'          sorted by position.
#' @param allNarrowPeaks a \code{GRanges} containing all narrow peaks from all 
#'          experiments sorted by position.
#' 
#' @return an object of \code{class} "commonFeatures". 
#' 
#' @author Astrid Louise Deschenes
#' @importFrom BiocGenerics start end
#' @importFrom stringr str_split
#' @importFrom IRanges IRanges 
#' @importFrom GenomicRanges GRanges findOverlaps seqinfo seqnames subjectHits
#' @keywords internal
findConsensusPeakRegionsForOneChrom <- function(chrName, extendingSize, 
                includeAllPeakRegion, minNbrExp, allPeaks, allNarrowPeaks) {
    
    # Subset peaks and narrow peaks using the specified chromosome name
    peaks <- sort(subset(allPeaks, as.character(seqnames(allPeaks)) == chrName))
    narrowPeaks <- sort(subset(allNarrowPeaks, 
                        as.character(seqnames(allNarrowPeaks)) == chrName))
    
    # Variable containing final consensus peak regions
    regions <- GRanges()
    
    if (length(peaks) > 0 && length(narrowPeaks) > 0) {
        # Variables initialization
        current <- NULL
        rightBoundary <- NULL
        namesVec <- vector()
        bad <- FALSE
        pos <- 1
        region_width <- 2 * extendingSize
        
        # All peak are tested
        # A primary region starting at peak position and of size 
        # 2 * extendingSize. All peaks included in the region are used to 
        # calcule the median of the peaks. Using the median position, a new
        # region of size 2 * extendingSize is created with the median at its
        # center. All peaks included in the region are used to 
        # calcule the median of the peaks. The iteration goes on as long as
        # the set of peaks is not stable and the inital peak is not part
        # of the region.
        # When a region is fixed. 
        repeat  {
            current <- peaks[pos]
            rightBoundaryNew <- start(current)
            seq_name <- as.character(seqnames(current))
            rightBoundary <- NULL
            set <- NULL
            setNew <- NULL
            noRegionFound <- FALSE
            repeat {
                set <- setNew
                rightBoundary <- rightBoundaryNew
                # Find peaks that overlaps the region
                overlaps <- findOverlaps(query = GRanges(seqnames = seq_name,
                                    ranges=c(IRanges(rightBoundary, 
                                    rightBoundary + region_width))),
                                    subject = peaks)
                setNew <- peaks[subjectHits(overlaps)]
                if (length(setNew) == 0 || !(current$name %in% setNew$name)) {
                    # The current peak is not included in the current region
                    # The region will not be selected
                    noRegionFound <- TRUE
                    break
                }
                
                # Use the median of the peaks to set the new right boundary
                rightBoundaryNew <- median(start(setNew)) - extendingSize
                # Stop loop when the overlaping peaks are stable or 
                # when no peaks is found
                if (!is.null(set) && (length(set) == length(setNew)) && 
                        all(set == setNew)) break
            }
            
            if (noRegionFound) {
                # Treat the next position
                pos <- pos + 1
            } else {
                # Keep region only when the number of different experiments 
                # present is reached
                short_names <- sapply(X = set$name, 
                                    function(x) stringr::str_split(string = x, 
                                            pattern = ".bam")[[1]][1])
                if (length(unique(short_names)) >= minNbrExp) {
                    # Create one final region using the narrow information 
                    # for each peak present
                    minPos <- rightBoundaryNew
                    maxPos <- minPos + region_width
                    if (includeAllPeakRegion) {
                        peakMedian <- rightBoundaryNew + extendingSize
                        for (name in unique(short_names)) {
                            peaksForOneExp <- set[short_names == name]
                        
                            closessPeaks <- which(abs(start(peaksForOneExp) - 
                                        peakMedian) == 
                                        min(abs(start(peaksForOneExp) - 
                                        peakMedian)))
                            peaksForOneExp <- peaksForOneExp[closessPeaks]
                            
                            newMax <- max(end(narrowPeaks[narrowPeaks$name 
                                                %in% peaksForOneExp$name]))
                            maxPos <- ifelse(newMax > maxPos, newMax, maxPos)
                        
                            newMin <- min(start(narrowPeaks[narrowPeaks$name 
                                                %in% peaksForOneExp$name]))
                            minPos <- ifelse(newMin < minPos, newMin, minPos)
                        }
                    }
                    
                    newRegion <- GRanges(seqnames = seq_name, 
                                            IRanges(minPos, maxPos))
                    regions <- append(regions, newRegion)
                    
                    # Update overlapping peaks
                    overlaps <- findOverlaps(query = newRegion, subject = peaks)
                    setNew <- peaks[subjectHits(overlaps)]
                    if (!(current$name %in% setNew$name)) {
                        # If the current peak is not included in the current 
                        # region, the region will not be selected
                        stop(paste0("The current treated peak should be in ", 
                                    "the selected region.\n"))
                    }
                    
                    # Treat the position following last peak present
                    # in new region
                    pos <- max(subjectHits(overlaps)) + 1
                } else {
                    pos <- pos + 1
                }
            }
        
            # Stop loop when all peaks are treated
            if (pos >= length(peaks)) break
        }
    }
    
    return(regions)
}
