###################################################
# Created by Astrid Louise Deschenes
# 2015-06-02
###################################################

###################################################
## Test the findConsensusPeakRegions.R functions
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "sharedBed" )
}

### }}}

data(Hosa_A549_FOSL2_01_Peaks)
data(Hosa_A549_FOSL2_01_NarrowPeaks)

###################################################
## Test the findConsensusPeakRegions() function
###################################################

## Test the result when a numerical is passed as narrowPeaks parameter
test.findConsensusPeakRegions_with_narrowPeaks_numerical <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 444, 
            peaks = GRanges()), error = conditionMessage)
    exp <- "narrowPeaks must be a GRanges object"
    message <- paste0("findConsensusPeakRegions_with_narrowPeaks_numerical() ",
            "- A numerical as narrowPeaks parameter did not generated ", 
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a integer is passed as peaks parameter
test.findConsensusPeakRegions_with_peaks_integer <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = GRanges(), 
            peaks = 444), error = conditionMessage)
    exp <- "peaks must be a GRanges object"
    message <- paste0("findConsensusPeakRegions_with_peaks_integer() ",
            "- A integer as peaks parameter did not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a list is passed as peaks parameter
test.findConsensusPeakRegions_with_peaks_list <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = GRanges(), 
            peaks = list()), error = conditionMessage)
    exp <- "peaks must be a GRanges object"
    message <- paste0("findConsensusPeakRegions_with_narrowPeaks_integer() ",
            "- A list as peaks parameter did not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a empty GRanges is passed as narrowPeaks parameter
test.findConsensusPeakRegions_with_narrowPeaks_empty_GRanges <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = GRanges(), 
            peaks = GRanges()), error = conditionMessage)
    exp <- "narrowPeaks must be a GRanges object with at least one entry"
    message <- paste0("findConsensusPeakRegions_with_narrowPeaks_empty_GRanges",
            "s() - A empty GRanges as narrowPeaks parameter did not ", 
            "generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when a empty GRanges is passed as peaks parameter
test.findConsensusPeakRegions_with_peaks_empty_GRanges <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks, peaks = GRanges()), 
            error = conditionMessage)
    exp <- "peaks must be a GRanges object with at least one entry"
    message <- paste0("findConsensusPeakRegions_with_peaks_empty_GRanges",
            "s() - A empty GRanges as peaks parameter did not generated ", 
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when two GRanges of different lengths are passed as
## parameters
test.findConsensusPeakRegions_with_diff_length_GRanges <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks[1:3], 
            peaks = Hosa_A549_FOSL2_01_Peaks[1:5]), 
            error = conditionMessage)
    exp <- "narrowPeaks and peaks must have the same number of elements"
    message <- paste0("findConsensusPeakRegions_with_diff_length_GRanges",
            "s() - Two GRanges of different lengths did not generated ", 
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when GRanges without names is passed as narrowPeaks
## parameter
test.findConsensusPeakRegions_narrowPeaks_with_no_name_GRanges <- function() {
    seqinfo <- Seqinfo(paste0("chr", 1:2), c(1000, 2000), NA, "mock1")
    gr <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(1, 1)),
                ranges = IRanges(1:2, width = 2:1, names=head(letters,2)),
                strand = Rle(strand(c("-", "+")), c(1, 1)),
                score = 1:2, GC = seq(1, 0, length=2),
                seqinfo=seqinfo)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = gr, 
                    peaks = Hosa_A549_FOSL2_01_Peaks[2:3]), 
                    error = conditionMessage)
    exp <- paste0("narrowPeaks and peaks must have defined names so that ", 
                  "each narrowPeaks entry can be associated to a peaks entry")
    message <- paste0("findConsensusPeakRegions_peaks_with_no_name_GRanges
                      ",
                      " - A GRanges without names used as narrowPeaks ", 
                      "parameter did not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when GRanges without names is passed as peaks
## parameter
test.findConsensusPeakRegions_peaks_with_no_name_GRanges <- function() {
    seqinfo <- Seqinfo(paste0("chr", 1:2), c(1000, 2000), NA, "mock1")
    gr <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(1, 1)),
                  ranges = IRanges(1:2, width = 2:1, names=head(letters,2)),
                  strand = Rle(strand(c("-", "+")), c(1, 1)),
                  score = 10:11, GC = seq(1, 0, length=2),
                  seqinfo=seqinfo)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                    Hosa_A549_FOSL2_01_Peaks[2:3], 
                    peaks = gr), 
                    error = conditionMessage)
    exp <- paste0("narrowPeaks and peaks must have defined names so that ", 
                  "each narrowPeaks entry can be associated to a peaks entry")
    message <- paste0("findConsensusPeakRegions_peaks_with_no_name_GRanges",
                      " - A GRanges without names used as peaks ", 
                      "parameter did not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when two GRanges with different names are passed as
## parameters
test.findConsensusPeakRegions_with_diff_names_GRanges <- function() {
    cl <- match.call()
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks[3:5], 
            peaks = Hosa_A549_FOSL2_01_Peaks[2:4]), 
            error = conditionMessage)
    exp <- paste0("All narrowPeaks entry must have an equivalent peaks ", 
                  "entry recognizable by a identical name")
    message <- paste0("findConsensusPeakRegions_with_diff_names_GRanges",
            " - Two GRanges with different names did not generated ", 
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

