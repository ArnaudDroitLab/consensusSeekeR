#' @title Parameter validation for the \code{\link{findCommonFeatures}} 
#'      function
#' 
#' @description Validation of all parameters needed by the public
#'      \code{\link{findCommonFeatures}} function.
#' 
#' @param narrowpeaksBEDFiles a \code{vector} containing the BED files to
#'          use for the regions selection.
#' @param chrList a \code{vector} containing the name of the chromosomes to 
#'          analyze or the name \code{"ALL"} which indicate that all
#'          chromosomes must be analyzed. When \code{NULL}, no
#'          new term is added. Default : \code{NULL}.
#' @param padding a \code{numeric}. Default = 250.
#' @param minNbrExp a \code{numeric} indicating the minimum number of BED files
#'          in which a peak must be present for a region to be retained. The
#'          numeric must be a positive value inferior or equal to the number of 
#'          files present in the \code{narrowpeaksBEDFiles} parameter.
#'          Default = 1.
#' @param nbrThreads a \code{numeric} indicating the number of threads to use
#'          in parallel.
#' 
#' @return \code{0} indicating that all parameters validations have been
#'      successful.
#' 
#' @author Astrid Louise Deschenes
#' @keywords internal
findCommonFeaturesValidation <- function(narrowpeaksBEDFiles, chrList, 
                                         padding, minNbrExp, nbrThreads) {
    
    if (is.vector(narrowpeaksBEDFiles) && (!is.character(narrowpeaksBEDFiles) 
            || !all(sapply(narrowpeaksBEDFiles, file.exists)))) {
        stop("peaksBEDlist must be a vector of existing BED files")
    }
    
    if (chrList != "ALL" && !(is.vector(chrList) && is.character(chrList))) {
        stop(paste0("chrList must be either be the value \"ALL\" or a ",
             "vector of chromosomes names"))
    }
      
    if (!isInteger(padding) || padding < 1 ) {
        stop("padding must be a non-negative integer")
    }
    
    if (!isInteger(minNbrExp) || minNbrExp < 1  || 
            minNbrExp > length(narrowpeaksBEDFiles)) {
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
#' @param chrName an object of \code{class} "formula" which contains a symbolic
#'          model formula.
#' @param padding a \code{data.frame} containing the variables in the model.
#' @param minNbrExp
#' @param allPeaks a \code{GRanges} indicating which column from the 
#'          \code{data} must be added to the formula. When \code{NULL}, no
#'          new term is added. Default : \code{NULL}.
#' @param allNarrowPeaks a \code{GRanges}.
#' 
#' @return an object of \code{class} "commonFeatures". 
#' 
#' @author Astrid Louise Deschenes
#' @importFrom BiocGenerics start end
#' @importFrom stringr str_split
#' @importFrom IRanges IRanges 
#' @importFrom GenomicRanges GRanges findOverlaps seqinfo seqnames subjectHits
#' @keywords internal
findCommonFeaturesForOneChrom <- function(chrName, padding, minNbrExp, 
                        allPeaks, allNarrowPeaks) {
    
    ## Only keep data related to the selected chromosome
    peaks <- sort(subset(allPeaks, seqnames(allPeaks) == chrName))
    narrowPeaks <- sort(subset(allNarrowPeaks, 
                            seqnames(allNarrowPeaks) == chrName))
    
    # Variables initialization
    regions <- GRanges()
    strangeRegions <- GRanges()
    current <- NULL
    rightBoundary <- NULL
    namesVec <- vector()
    bad <- FALSE
    pos <- 1
    
    print(paste("Length(peaks):", length(peaks)))
    
    repeat  {
        current <- peaks[pos]
        rightBoundaryNew <- start(current)
        seq_name <- as.character(seqnames(current))
        rightBoundary <- NULL
        set <- NULL
        setNew <- NULL
        bad <- FALSE
        region_width <- 2 * padding
        repeat {
            set <- setNew
            rightBoundary <- rightBoundaryNew
            # Find peaks that overlaps the region
            overlaps <- findOverlaps(query = GRanges(seqnames = seq_name,
                            ranges=c(IRanges(rightBoundary, 
                            rightBoundary + region_width))),
                            subject = peaks)
            setNew <- peaks[subjectHits(overlaps)]
            if (!(current$name %in% setNew$name)) {
                # The current peak is not included in the current region
                # The region will not be selected
                #strangeRegions <- append(strangeRegions, current)
                bad <- TRUE
                break
            }
            
            # Use the median of the peaks to set the new right boundary
            rightBoundaryNew <- median(start(setNew)) - padding
            # Stop loop when the overlaping peaks are stable or 
            # when no peaks are found
            if (!is.null(set) && (length(set) == length(setNew)) && 
                    all(set == setNew)) break
        }
        
        if (bad) {
            # Treat the next position
            pos <- pos + 1
        } else {
            # Keep region only when the number of different experiments present
            # is reached
            short_names <- sapply(X = set$name, 
                            function(x) stringr::str_split(string = x, 
                                            pattern = ".bam")[[1]][1])
            if (length(unique(short_names)) > minNbrExp) {
                # Create one final region using the narrow information 
                # for each peak present
                minPos <- rightBoundaryNew
                peakMedian <- rightBoundaryNew + padding
                maxPos <- peakMedian + padding
                for (i in unique(short_names)) {
                    peaksForOneExp <- set[short_names == i]
                    
                    closessPeak <- which(abs(start(peaksForOneExp) - peakMedian)
                                            == 
                                    min(abs(start(peaksForOneExp)- peakMedian)))
                
                    firstPeak <- peaksForOneExp[closessPeak[1]]
                    newMax <- end(narrowPeaks[narrowPeaks$name %in% 
                                                    lastPeak$name])
                    maxPos <- ifelse(newMax > maxPos, newMax, maxPos)
                    
                    lastPeak <- peaksForOneExp[closessPeak[length(closessPeak)]]
                    newMin <- start(narrowPeaks[narrowPeaks$name %in%  
                                                firstPeak$name])
                    minPos <- ifelse(newMin < minPos, newMin, minPos)
                }
                
                newRegion <- GRanges(seqnames = seq_name, 
                                        IRanges(minPos, maxPos))
                regions <- append(regions, newRegion)
                
                # Update overlapping peaks
                overlaps <- findOverlaps(query = newRegion, subject = peaks)
                setNew <- peaks[subjectHits(overlaps)]
                if (!(current$name %in% setNew$name)) {
                    # The current peak is not included in the current region
                    # The region will not be selected
                    stop(paste0("The current treated peak should be in the ", 
                                    "selected region.\n"))
                }
                
                # Treat the position following last peak present in new region
                pos <- max(subjectHits(overlaps)) + 1
            } else {
                pos <- pos + 1
            }
        }
        # Stop loop when all peaks are treated
        if (pos >= length(peaks)) break
    }
    result <- list(features=regions)
    class(result) <- "commonFeatures"
    return()
}
