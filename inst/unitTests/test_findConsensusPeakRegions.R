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
data(Hosa_A549_FOSL2_01_Peaks_partial)
data(Hosa_A549_FOSL2_01_NarrowPeaks_partial)
data(Hosa_A549_FOXA1_01_Peaks)
data(Hosa_A549_FOXA1_01_NarrowPeaks)
data(Hosa_A549_FOXA1_01_Peaks_partial)
data(Hosa_A549_FOXA1_01_NarrowPeaks_partial)

###########################################################
## Test the findConsensusPeakRegions() function parameters
###########################################################

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
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:3], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:5]), 
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
                    peaks = Hosa_A549_FOSL2_01_Peaks_partial[2:3]), 
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
                    Hosa_A549_FOSL2_01_Peaks_partial[2:3], 
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
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[3:5], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[2:4]), 
            error = conditionMessage)
    exp <- paste0("All narrowPeaks entry must have an equivalent peaks ", 
                  "entry recognizable by a identical name")
    message <- paste0("findConsensusPeakRegions_with_diff_names_GRanges",
            " - Two GRanges with different names did not generated ", 
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrList is a numerical
test.findConsensusPeakRegions_with_numerical_chrList <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], chrList = 444),
            error = conditionMessage)
    exp <- paste0("chrList must either be the value \"ALL\" or a ",
                  "vector of chromosomes names")
    message <- paste0("findConsensusPeakRegions_with_strange_chrList",
                      " - Numerical as chrList did not generated ", 
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrList is strange string
test.findConsensusPeakRegions_with_strange_chrList <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], chrList = "ALLO"),
            error = conditionMessage)
    exp <- "No chromosome correspond to the given parameter: ALLO"
    message <- paste0("findConsensusPeakRegions_with_strange_chrList",
                " - Strange string as chrList did not generated ", 
                "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrList is a list with not corresponding name
test.findConsensusPeakRegions_with_list_strange_name_as_chrList <- function() {
    testList <- c("ALLO", "BYE")
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], chrList = testList),
            error = conditionMessage)
    exp <- paste0("No chromosome correspond to the given parameters: ", 
                  paste0(testList,  collapse = ", "))
    message <- paste0("findConsensusPeakRegions_with_list_strange_name_as",
            "_chrList - List of strange string as chrList did not generated ", 
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when zero as extendingSize
test.findConsensusPeakRegions_with_zero_as_extendingSize <- function() {
    testList <- c("ALLO", "BYE")
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], extendingSize = 0),
            error = conditionMessage)
    exp <- "extendingSize must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_with_zero_as_extendingSize",
                      " - Xero as extendingSize did not generated ", 
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}
    
## Test the result when negative as extendingSize
test.findConsensusPeakRegions_with_negative_as_extendingSize <- function() {
    testList <- c("ALLO", "BYE")
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], extendingSize = -90),
            error = conditionMessage)
    exp <- "extendingSize must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_with_negative_as_extendingSize",
                    " - Negative as extendingSize did not generated ", 
                    "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when string as extendingSize
test.findConsensusPeakRegions_with_string_as_extendingSize <- function() {
    testList <- c("ALLO", "BYE")
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
            extendingSize = "444"), error = conditionMessage)
    exp <- "extendingSize must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_with_string_as_extendingSize",
                      " - String as extendingSize did not generated ", 
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when string as includeAllPeakRegion
test.findConsensusPeakRegions_string_as_includeAllPeakRegion <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
            includeAllPeakRegion = "444"), error = conditionMessage)
    exp <- "includeAllPeakRegion must be a logical value"
    message <- paste0("findConsensusPeakRegions_string_as_includeAllPeakRegion",
                " - String as includeAllPeakRegion did ", 
                "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when numerical as includeAllPeakRegion
test.findConsensusPeakRegions_numerical_as_includeAllPeakRegion <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            c(Hosa_A549_FOSL2_01_NarrowPeaks_partial, 
            Hosa_A549_FOXA1_01_NarrowPeaks_partial), 
            peaks = c(Hosa_A549_FOSL2_01_Peaks_partial, 
            Hosa_A549_FOXA1_01_Peaks_partial), 
            includeAllPeakRegion=333), error = conditionMessage)
    exp <- "includeAllPeakRegion must be a logical value"
    message <- paste0("findConsensusPeakRegions_numerical_as_",
                      "includeAllPeakRegion - Numerical as ", 
                      "includeAllPeakRegion did not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when string as minNbrExp
test.findConsensusPeakRegions_string_as_minNbrExp <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                minNbrExp = "444"), error = conditionMessage)
    exp <- "minNbrExp must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_string_as_minNbrExp",
                      " - String as minNbrExp did ", 
                      "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when zero as minNbrExp
test.findConsensusPeakRegions_zero_as_minNbrExp <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                minNbrExp = 0), error = conditionMessage)
    exp <- "minNbrExp must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_zero_as_minNbrExp",
                      " - Zero as minNbrExp did ", 
                      "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when negative as minNbrExp
test.findConsensusPeakRegions_negative_as_minNbrExp <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                                minNbrExp = -1), error = conditionMessage)
    exp <- "minNbrExp must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_negative_as_minNbrExp",
                      " - Negative as minNbrExp did ", 
                      "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when numerical as minNbrExp
test.findConsensusPeakRegions_numerical_as_minNbrExp <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                                minNbrExp = 9.3), error = conditionMessage)
    exp <- "minNbrExp must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_numerical_as_minNbrExp",
                      " - Numerical as minNbrExp did ", 
                      "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when list of integers as minNbrExp
test.findConsensusPeakRegions_list_of_integers_as_minNbrExp <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                              Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                              peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                              minNbrExp= c(9L, 3L)), error = conditionMessage)
    exp <- "minNbrExp must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_list_of_integers_as_minNbrExp",
                      " - List of integers as minNbrExp did ", 
                      "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}


## Test the result when zero as nbrThreads
test.findConsensusPeakRegions_zero_as_nbrThreads<- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                                nbrThreads = 0), error = conditionMessage)
    exp <- "nbrThreads must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_zero_as_nbrThreads",
                      " - Zero as nbrThreads did ", 
                      "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when negative as nbrThreads
test.findConsensusPeakRegions_negative_as_nbrThreads <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                                nbrThreads = -1), error = conditionMessage)
    exp <- "nbrThreads must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_negative_as_nbrThreads",
                      " - Negative as nbrThreads did ", 
                      "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when numerical as nbrThreads
test.findConsensusPeakRegions_numerical_as_nbrThreads <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                                nbrThreads = 9.3), error = conditionMessage)
    exp <- "nbrThreads must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_numerical_as_nbrThreads",
                      " - Numerical as nbrThreads did ", 
                      "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when list of integers as nbrThreads
test.findConsensusPeakRegions_list_of_integers_as_nbrThreads <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                            nbrThreads= c(9L, 3L)), error = conditionMessage)
    exp <- "nbrThreads must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_list_of_integers_as_nbrThreads",
                      " - List of integers as nbrThreads did ", 
                      "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}


###########################################################
## Test the findConsensusPeakRegions() function results
###########################################################

## Test the result when only one chromosome in list
test.findConsensusPeakRegions_for_one_chromosome <- function() {
    seqinfo <- Seqinfo(paste0("chr", 1), NA, NA, NA)
    exp <- GRanges(seqnames = Rle(c("chr1"), c(1)),
                  ranges = IRanges(start = c(249119914, 249120334, 249123074, 
                                             249132040, 249133011, 249152098, 
                                             249152823, 249153205, 249157198,
                                             249167214, 249167809, 249199968),
                                   end = c(249120424, 249121174, 249123574, 
                                           249132673, 249133517, 249152644, 
                                           249153397, 249153705, 249157698, 
                                           249167714, 249168802, 249200468)),
                  seqinfo=seqinfo)
    obs <- findConsensusPeakRegions(narrowPeaks = 
                c(Hosa_A549_FOSL2_01_NarrowPeaks_partial, 
                Hosa_A549_FOXA1_01_NarrowPeaks_partial), 
                peaks = c(Hosa_A549_FOSL2_01_Peaks_partial, 
                Hosa_A549_FOXA1_01_Peaks_partial), chrList="chr1")
    checkEquals(obs$consensusRanges, exp, msg = message)
}

