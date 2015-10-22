#' @title Parameter validation for the \code{\link{findConsensusPeakRegions}}
#' function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{findConsensusPeakRegions}} function.
#'
#' @param narrowPeaks a \code{GRanges} representing called peaks of signal
#' for all experiments.
#'
#' @param peaks a \code{GRanges} representing peaks for all experiments.
#'
#' @param chrList a \code{Seqinfo} containing the name and the length of the
#' chromosomes to analyze which indicate that all chromosomes must
#' be analyzed.
#'
#' @param extendingSize a \code{numeric} value indicating the size of padding
#' at each side of the peaks median position to create the consensus
#' region. The minimum size of the consensus region will be equal to
#' twice the value of the \code{extendingSize} parameter. The size of
#' the \code{extendingSize} must be a positive integer. Default = 250.
#'
#' @param expandToFitPeakRegion a \code{logical} indicating if the region set
#' by the \code{extendingSize} parameter is extended to include all
#' region of the peak closest to the peaks median position for each
#' experiment.
#'
#' @param shrinkToFitPeakRegion a \code{logical} indicating if the region size,
#' which is set by the \code{extendingSize} parameter is shrinked to
#' fit the narrow peak regions of the peaks when all those regions
#' are smaller than the consensus region.
#'
#' @param minNbrExp a \code{numeric} indicating the minimum number of BED files
#' in which a peak must be present for a region to be retained. The numeric
#' must be a positive value inferior or equal to the number of files
#' present in the \code{narrowPeakFiles} parameter. Default = 1.
#'
#' @param nbrThreads a \code{numeric} indicating the number of threads to use
#' in parallel.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
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
#'                         length(A549_FOXA1_01_Peaks_partial))
#' names(A549_FOXA1_01_NarrowPeaks_partial) <- rep("FOXA1_01",
#'                         length(A549_FOXA1_01_NarrowPeaks_partial))
#' names(A549_FOSL2_01_Peaks_partial) <- rep("FOSL2_01",
#'                         length(A549_FOSL2_01_Peaks_partial))
#' names(A549_FOSL2_01_NarrowPeaks_partial) <- rep("FOSL2_01",
#'                         length(A549_FOSL2_01_NarrowPeaks_partial))
#'
#' chrList <- Seqinfo("chr10", 135534747, NA)
#'
#' consensusSeekeR:::findConsensusPeakRegionsValidation(
#'     narrowPeaks = c(A549_FOXA1_01_NarrowPeaks_partial,
#'             A549_FOSL2_01_NarrowPeaks_partial),
#'     peaks = c(A549_FOXA1_01_Peaks_partial,
#'             A549_FOSL2_01_Peaks_partial),
#'     chrList = chrList,
#'     extendingSize = 110,
#'     expandToFitPeakRegion = FALSE,
#'     shrinkToFitPeakRegion = TRUE,
#'     minNbrExp = 2,
#'     nbrThreads = 1)
#'
#' @author Astrid Deschenes
#' @importFrom GenomeInfoDb Seqinfo seqinfo seqlengths
#' @keywords internal
findConsensusPeakRegionsValidation <- function(narrowPeaks, peaks, chrList,
            extendingSize, expandToFitPeakRegion, shrinkToFitPeakRegion,
            minNbrExp, nbrThreads) {

    if (!is.logical(expandToFitPeakRegion)) {
        stop("expandToFitPeakRegion must be a logical value")
    }

    if (!is.logical(shrinkToFitPeakRegion)) {
        stop("shrinkToFitPeakRegion must be a logical value")
    }

    if (!is(peaks, "GRanges")) {
        stop("peaks must be a GRanges object")
    }

    if (length(peaks) < 1) {
        stop("peaks must be a GRanges object with at least one entry")
    }

    if (!is(chrList, "Seqinfo")) {
        stop("chrList must be a Seqinfo object")
    }

    if (expandToFitPeakRegion | shrinkToFitPeakRegion) {
        ## narrowPeaks is only validated when used by code
        if (!is(narrowPeaks, "GRanges")) {
            stop("narrowPeaks must be a GRanges object")
        }

        if (length(narrowPeaks) < 1) {
            stop("narrowPeaks must be a GRanges object with at least one entry")
        }

        if(length(narrowPeaks) != length(peaks)) {
            stop("narrowPeaks and peaks must have the same number of elements")
        }

        if(is.null(narrowPeaks$name) || is.null(peaks$name)) {
            stop(paste0("narrowPeaks and peaks must have defined metadata name ",
                        "so that each narrowPeaks entry can be associated to ",
                        "a peaks entry"))
        }

        if(is.null(names(narrowPeaks)) || is.null(names(peaks))) {
            stop(paste0("narrowPeaks and peaks must have defined row names ",
                        "so that each entry can be associated to an experiment"))
        }

        if (!all(sort(narrowPeaks$name) == sort(peaks$name)) ||
                !all(sort(names(narrowPeaks)) == sort(names(peaks)))) {
            stop(paste0("All narrowPeaks entry must have an equivalent peaks ",
                        "entry recognizable by both an identical metadata name and an ",
                        "identical row name"))
        }
    }

    if (any(is.na(seqlengths(chrList)))) {
        stop("At least one chromosome length is missing in chrList")
    }

    if (!all(names(chrList) %in% names(seqinfo(peaks)))) {
        not_present <- names(chrList)[!(names(chrList)
                                        %in% names(seqinfo(peaks)))]
        if (length(not_present) == length(names(chrList))) {
            stop("No chromosome name from chrList is present in peak")
        }
    }

    if (!isInteger(extendingSize) || extendingSize < 1 ) {
        stop("extendingSize must be a non-negative integer")
    }

    if (!isInteger(minNbrExp) || minNbrExp < 1) {
        stop("minNbrExp must be a non-negative integer")
    }

    if (minNbrExp > length(unique(names(peaks)))) {
        stop(paste0("minNbrExp must be inferior or equal to the number ",
                    "of experiments presents in peaks. The number of",
                    "experiments is known by the number of differents row names ",
                    "in peaks."))
    }

    if (!isInteger(nbrThreads) || nbrThreads < 1 ) {
        stop("nbrThreads must be a non-negative integer")
    }

    return(0)
}

#' @title Validate if a value is an integer
#'
#' @description Validate if the value passed to the function is an integer or
#' not. To be considered as an integer, the value must have a length
#' of 1. The type of value can be a \code{integer} or \code{numerical}.
#' However, a \code{numerical} must have the same value
#' once casted to a \code{integer}.  A \code{vector} of
#' integers will returned \code{FALSE}.
#'
#' @param value an object to validate.
#'
#' @return \code{TRUE} is the parameter is a integer; otherwise \code{FALSE}
#'
#' @author Astrid Deschenes
#' @keywords internal
isInteger <- function(value) {
    return((is.integer(value) && length(value) == 1) || (is.numeric(value) &&
                as.integer(value) == value) && length(value) == 1)
}

#' @title Extract regions sharing features in more than one experiment for
#' one specific chromosome.
#'
#' @description Find regions sharing the same features for a minimum number of
#' experiments using called peaks of signal enrichment based on
#' pooled, normalized data (mainly coming from narrowPeak files). Tje
#' analysis is limited to one chromosome.
#' The peaks and narrow peaks are used to identify
#' the consensus regions. The minimum number of experiments that must
#' have at least on peak in a region so that it is retained as a
#' consensus region is specified by user, as well as the size of
#' mining regions.
#'
#' @param chrName a \code{character}, the name of the chromosome to analyze.
#'
#' @param allPeaks a \code{GRanges} containing all peaks from all experiments
#' sorted by position.
#'
#' @param allNarrowPeaks a \code{GRanges} containing all narrow peaks from all
#' experiments sorted by position.
#'
#' @param extendingSize a \code{numeric} value indicating the size of padding
#' at each side of the peaks median position to create the consensus
#' region. The minimum size of the consensu region will be equal to
#' twice the value of the \code{extendingSize} parameter. The size of
#' the \code{extendingSize} must be a positive integer. Default = 250.
#'
#' @param expandToFitPeakRegion a \code{logical} indicating if the region set by
#' the \code{extendingSize} parameter is extended to include all narrow peak
#' regions. Only the narrow peaks regions of the peaks included in the
#' unextended region are used during the extension process. It is possible
#' that has a side effect, adding peaks are being included in the
#' extended region.
#'
#' @param shrinkToFitPeakRegion a \code{logical} indicating if the region set
#' by the \code{extendingSize} parameter is shrinked to fit the narrow
#' peak regions.
#'
#' @param minNbrExp a \code{numeric} indicating the minimum number of BED files
#' in which a peak must be present for a region to be retained. The
#' numeric must be a positive value inferior or equal to the number
#' of files present in the \code{narrowPeakFiles} parameter. Default = 1.
#'
#' @param chrList a \code{Seqinfo} containing the name and the length of the
#' chromosomes to analyze.
#'
#' @return an object of \code{class} "commonFeatures".
#'
#' @author Astrid Deschenes
#' @importFrom BiocGenerics start end
#' @importFrom stringr str_split
#' @importFrom IRanges IRanges median ranges "ranges<-"
#' @importFrom GenomicRanges GRanges findOverlaps seqnames subjectHits ranges
#' @importFrom GenomeInfoDb Seqinfo
#' @keywords internal
findConsensusPeakRegionsForOneChrom <- function(chrName, allPeaks,
                allNarrowPeaks, extendingSize, expandToFitPeakRegion,
                shrinkToFitPeakRegion,
                minNbrExp, chrList) {

    # Detect if narrowPeaks are needed or not
    areNarrowPeaksUsed <- expandToFitPeakRegion | shrinkToFitPeakRegion

    # Subset peaks using the specified chromosome name
    peaks <- sort(subset(allPeaks,
                        as.character(seqnames(allPeaks)) == chrName))

    # Subset narrow peaks using the spcified chromosome name when necessary
    if (areNarrowPeaksUsed) {
        narrowPeaks <- sort(subset(allNarrowPeaks,
                        as.character(seqnames(allNarrowPeaks)) == chrName))
    } else {
        narrowPeaks <- GRanges()
    }

    # Get info for the specific chromosome
    chrInfo <- chrList[chrName]

    # GRanges containing final consensus peak regions
    maxLength <- 10000
    increment <- 10000
    regions   <- GRanges(seqnames = Rle(chrName, maxLength),
                            rep(IRanges(1, 1), maxLength))


    nbrPeaks   <- length(peaks)
    nbrRegions <- 0L

    if (nbrPeaks > 0) {
        # Variables initialization
        current <- NULL
        rightBoundary <- NULL
        bad <- FALSE
        pos <- 1
        region_width <- 2 * extendingSize

        tempGRange <- GRanges(seqnames = chrName, IRanges(1, 1))

        # All peak are tested.
        # A primary region starting at peak position and of size
        # 2 * extendingSize. All peaks included in the region are used to
        # calcule the median of the peaks. Using the median position, a new
        # region of size 2 * extendingSize is created with the median at its
        # center. All peaks included in the region are used to
        # calcule the median of the peaks. The iteration goes on as long as
        # the set of peaks is not stable and the inital peak is not part
        # of the region.
        # When a region is fixed, use narrowPeak information to ajust
        # boundaries if needed.
        repeat  {
            current <- peaks[pos]
            rightBoundaryNew <- start(current)
            rightBoundary <- NULL
            set <- NULL
            setNew <- NULL
            noRegionFound <- FALSE
            repeat {
                set <- setNew
                rightBoundary <- rightBoundaryNew

                # Update GRange to fit the new region
                ranges(tempGRange) <- IRanges(rightBoundary,
                                             rightBoundary + region_width)

                # Find peaks that overlap the new region
                overlaps <- findOverlaps(query = tempGRange, subject = peaks)
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
                if (length(unique(names(set))) >= minNbrExp) {
                    # Create one final region using the narrow information
                    # for each peak present
                    minPos <- rightBoundaryNew
                    maxPos <- minPos + region_width

                    if (areNarrowPeaksUsed) {
                        narrowPeaksSet <- narrowPeaks[narrowPeaks$name %in%
                                                          set$name]

                        minLeft <- min(start(narrowPeaksSet))
                        maxRight <- max(end(narrowPeaksSet))

                        if (expandToFitPeakRegion) {
                            minPos <- ifelse(minLeft < minPos, minLeft, minPos)
                            maxPos <- ifelse(maxRight > maxPos,
                                                maxRight, maxPos)
                        }

                        if (shrinkToFitPeakRegion) {
                            minPos <- ifelse(minLeft > minPos, minLeft, minPos)
                            maxPos <- ifelse(maxRight < maxPos,
                                                maxRight, maxPos)
                        }
                    }

                    # Validate that minimum position is not negative
                    if (minPos < 1) {
                        minPos <- 1
                    }

                    # Validate that maximum position is not
                    # superior to chromosome length
                    if (maxPos > seqlengths(chrInfo)) {
                        maxPos <- seqlengths(chrInfo)
                    }

                    # Update total number of consensus regions
                    nbrRegions <- nbrRegions + 1L

                    # Adapt size of GRanges containing final results
                    # when too many regions
                    if (nbrRegions > maxLength) {
                        regions   <- append(regions, GRanges(seqnames =
                                            Rle(chrName, increment),
                                            rep(IRanges(1, 1), increment)))
                        maxLength <- maxLength + increment
                    }

                    # Update GRanges to contain the value of the new region
                    ranges(regions[nbrRegions])<- IRanges(minPos, maxPos)

                    # Update overlapping peaks
                    overlaps <- findOverlaps(query = regions[nbrRegions],
                                                subject = peaks)

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
            if (pos >= nbrPeaks) break
        }
    }

    # Adjust size of returned GRanges to the real number of consensus regions
    if (nbrRegions == 0) {
        regions <- GRanges()
    } else {
        regions <- regions[1:nbrRegions]
    }

    return(regions)
}
