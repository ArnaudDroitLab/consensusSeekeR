#' @title TODO
#' 
#' @description TODO
#' 
#' @param chrName an object of \code{class} "formula" which contains a symbolic
#'          model formula.
#' @param padding a \code{data.frame} containing the variables in the model.
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
#' @importFrom GenomicRanges GRanges
#' @keywords internal
findCommonFeaturesValidation <- function(peaksBEDlist, narrowpeaksBEDlist, 
                                    chrList, padding = 500, minNbrExp = 1) {
    
    if (is.vector(peaksBEDlist) && (!is.character(peaksBEDlist) || 
                !all(sapply(peaksBEDlist, file.exists)))) {
        stop("peaksBEDlist must be a vector of existing BED files")
    }
    
    if (is.vector(narrowpeaksBEDlist) && 
                (!is.character(narrowpeaksBEDlist) || 
                !all(sapply(peaksBEDlist, file.exists)))) {
            stop("peaksBEDlist must be a vector of existing BED files")
    }
        
    if (!isInteger(padding) || padding < 1 ) {
        stop("padding must be a non-negative integer")
    }
    
    if (!isInteger(minNbrExp) || minNbrExp < 1  || 
            minNbrExp > length(peakBEDlist)) {
        stop(paste0("minNbrExp must be a non-negative integer inferior or ", 
                "equal to the number of experiments."))
    }
      
}

#' @title Validate if a value is an integer
#' 
#' @description TODO
#'
#' @param value an object to validate.
#' 
#' @return 
#' 
#' @author Astrid Louise Deschenes
#' @keywords internal
isInteger <- function(value) {
    return((is.numeric(minNbrExp) || is.integer(cores)) && 
                as.integer(cores) == cores)
}

#' @title TODO
#' 
#' @description TODO
#' 
#' @param chrName an object of \code{class} "formula" which contains a symbolic
#'          model formula.
#' @param padding a \code{data.frame} containing the variables in the model.
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
#' @importFrom GenomicRanges GRanges
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
    
    repeat  {
        current <- peaks[pos]
        rightBoundaryNew <- start(current)
        rightBoundary <- NULL
        set <- NULL
        setNew <- NULL
        bad <- FALSE
        repeat {
            set <- setNew
            rightBoundary <- rightBoundaryNew
            # Find peaks that overlaps the region
            overlaps <- findOverlaps(query = GRanges(seqnames = c(as.character(seqnames(current))),
                            ranges=c(IRanges(rightBoundary, rightBoundary+(2*padding)))),
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
            rightBoundaryNew <- median(start(setNew))-padding
            # Stop loop when the overlaping peaks are stable or when no peaks are found
            if (!is.null(set) && (length(set) == length(setNew)) && all(set == setNew)) break
        }
        
        if (bad) {
            # Treat the next position
            pos <- pos + 1
        } else {
            # Keep region only when peaks from more than one experience are present
            short_names <- sapply(X = set$name, function(x) stringr::str_split(string = x, pattern = ".bam")[[1]][1])
            if (length(unique(short_names)) > 1) {
                # Create one final region using the narrow information for each peak present
                minPos <- ifelse(min(BiocGenerics::start(narrowPeaks[narrowPeaks$name %in%  set$name])), rightBoundary)
                maxPos <- rightBoundary + padding
                for (i in unique(short_names)) {
                    peaksForOneExp <- set[short_names == i]
                    
                    firstPeak <- peaksForOneExp[1]
                    lastPeak <- peaksForOneExp[length(peaksForOneExp)]
                    
                    newMax <- BiocGenerics::end(narrowPeaks[narrowPeaks$name %in%  firstPeak$name])
                    maxPos <- ifelse(newMax > maxPos, newMax, maxPos)
                    
                    newMin <- BiocGenerics::start(narrowPeaks[narrowPeaks$name %in%  lastPeak$name])
                    minPos <- ifelse(newMin < minPos, newMin, minPos)
                }
                
                newRegion <- GRanges(seqnames = as.character(seqnames(set[1])), IRanges(minPos, maxPos))
                regions <- append(regions, newRegion)
                # Update overlapping peaks
                overlaps <- findOverlaps(query = newRegion, subject = peaks)
            }
            # Treat the position following last peak present in new region
            pos<-max(subjectHits(overlaps))+1
        }
        # Stop loop when all peaks are treated
        if (pos >= length(peaks)) break
    }
    result <- list(features=regions)
    class(result) <- "commonFeatures"
    return()
}
