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
    library( "consensusSeekeR" )
}

### }}}

data(Hosa_A549_FOSL2_01_Peaks_partial)
data(Hosa_A549_FOSL2_01_NarrowPeaks_partial)
data(Hosa_A549_FOXA1_01_Peaks_partial)
data(Hosa_A549_FOXA1_01_NarrowPeaks_partial)

names(Hosa_A549_FOXA1_01_Peaks_partial) <- 
    rep("Hosa_A549_FOXA1_01", length(Hosa_A549_FOXA1_01_Peaks_partial))
names(Hosa_A549_FOXA1_01_NarrowPeaks_partial) <-
    rep("Hosa_A549_FOXA1_01", length(Hosa_A549_FOXA1_01_NarrowPeaks_partial))
names(Hosa_A549_FOSL2_01_Peaks_partial) <-
    rep("Hosa_A549_FOSL2_01", length(Hosa_A549_FOSL2_01_Peaks_partial))
names(Hosa_A549_FOSL2_01_NarrowPeaks_partial) <-
    rep("Hosa_A549_FOSL2_01", length(Hosa_A549_FOSL2_01_NarrowPeaks_partial))

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
            Hosa_A549_FOSL2_01_NarrowPeaks_partial, peaks = GRanges()), 
            error = conditionMessage)
    exp <- "peaks must be a GRanges object with at least one entry"
    message <- paste0(" findConsensusPeakRegions_with_peaks_empty_GRanges",
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
    message <- paste0(" findConsensusPeakRegions_with_diff_length_GRanges",
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
    exp <- paste0("narrowPeaks and peaks must have defined metadata name so ", 
                  "that each narrowPeaks entry can be associated to a ", 
                  "peaks entry")
    message <- paste0("findConsensusPeakRegions_peaks_with_no_name_GRanges",
                      " - A GRanges without names used as narrowPeaks ", 
                      "parameter did not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when GRanges without row name is passed as narrowPeaks
## parameter
test.findConsensusPeakRegions_narrowPeaks_with_no_row_name_GRanges <- function() {
    seqinfo <- Seqinfo(paste0("chr", 1:2), c(1000, 2000), NA, "mock1")
    gr <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(1, 1)),
                ranges = IRanges(1:2, width = 2:1),
                strand = Rle(strand(c("-", "+")), c(1, 1)),
                score = 1:2, GC = seq(1, 0, length=2), seqinfo=seqinfo)
    gr$name = paste0("peak", 1:2)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = gr, 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[2:3]), 
            error = conditionMessage)
    exp <- paste0("narrowPeaks and peaks must have defined row names ", 
                  "so that each entry can be associated to an ", 
                  "experiment")
    message <- paste0(" findConsensusPeakRegions_narrowPeaks_with_no_row_name_GRanges",
                      " - A GRanges without row name used as narrowPeaks ", 
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
                  score = 10:11, GC = seq(1, 0, length=2), seqinfo=seqinfo)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                    Hosa_A549_FOSL2_01_Peaks_partial[2:3], 
                    peaks = gr), error = conditionMessage)
    exp <- paste0("narrowPeaks and peaks must have defined metadata name ", 
            "so that each narrowPeaks entry can be associated to ", 
            "a peaks entry")
    message <- paste0(" findConsensusPeakRegions_peaks_with_no_name_GRanges",
                      " - A GRanges without names used as peaks ", 
                      "parameter did not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when GRanges without row name is passed as peaks
## parameter
test.findConsensusPeakRegions_peaks_with_no_row_name_GRanges <- function() {
    seqinfo <- Seqinfo(paste0("chr", 1:2), c(1000, 2000), NA, "mock1")
    gr <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(1, 1)),
                  ranges = IRanges(1:2, width = 2:1),
                  strand = Rle(strand(c("-", "+")), c(1, 1)),
                  score = 10:11, GC = seq(1, 0, length=2), seqinfo=seqinfo)
    gr$name <- paste0("peak", 1:2)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                            Hosa_A549_FOSL2_01_Peaks_partial[2:3], 
                            peaks = gr), error = conditionMessage)
    exp <- paste0("narrowPeaks and peaks must have defined row names ", 
                  "so that each entry can be associated to ", 
                  "an experiment")
    message <- paste0(" findConsensusPeakRegions_peaks_with_no_row_name_",
                      "GRanges - A GRanges without names used as peaks ", 
                      "parameter did not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when two GRanges with different names are passed as
## parameters
test.findConsensusPeakRegions_with_diff_names_GRanges <- function() {
    seqinfo <- Seqinfo(paste0("chr", 1:2), c(1000, 2000), NA, "mock1")
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[3:5], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[2:4]), chrInfo = seqinfo, 
            error = conditionMessage)
    exp <- paste0("All narrowPeaks entry must have an equivalent peaks ", 
        "entry recognizable by both an identical metadata name and an ", 
        "identical row name")
    message <- paste0("findConsensusPeakRegions_with_diff_names_GRanges",
            " - Two GRanges with different names did not generated ", 
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrList is a numerical
test.findConsensusPeakRegions_with_numerical_chrList <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], chrInfo = 444),
            error = conditionMessage)
    exp <- paste0("chrList must be a Seqinfo object")
    message <- paste0(" findConsensusPeakRegions_with_strange_chrList",
                    " - Numerical as chrList did not generated ", 
                    "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrList is strange string
test.findConsensusPeakRegions_with_strange_chrList <- function() {
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], chrInfo = "ALLO"),
            error = conditionMessage)
    exp <- "chrList must be a Seqinfo object"
    message <- paste0("findConsensusPeakRegions_with_strange_chrList",
                " - Strange string as chrList did not generated ", 
                "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrList is a list 
test.findConsensusPeakRegions_with_list_strange_name_as_chrList <- function() {
    testList <- c("ALLO", "BYE")
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], chrInfo = testList),
            error = conditionMessage)
    exp <- "chrList must be a Seqinfo object"
    message <- paste0("findConsensusPeakRegions_with_list_strange_name_as",
            "_chrList - List of strange string as chrList did not generated ", 
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrList with a missing length
test.findConsensusPeakRegions_with_missing_length_in_chrList <- function() {
    chrList <- Seqinfo(paste0("chr", c(1,2)), c(NA, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                    Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                    peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                    chrInfo = chrList),
                    error = conditionMessage)
    exp <- paste0("At least one chromosome length is missing in chrList")
    message <- paste0(" findConsensusPeakRegions_with_missing_length_",
                      "in_chrList - Absent chromosome in chrList did ", 
                      "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrList with a non existing chromosome
test.findConsensusPeakRegions_with_absent_chr_as_chrList <- function() {
    chrList <- Seqinfo(paste0("chr", c(1,40)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                        Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                        peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                        chrInfo = chrList),
                        error = conditionMessage)
    exp <- paste0("At least one chromosome name present in chrList is ",
            "not present in peak : chr40")
    message <- paste0(" findConsensusPeakRegions_with_absent_chr_as_chrList",
            " - Absent chromosome in chrList did not generated ", 
            "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrList with a non existing chromosomes
test.findConsensusPeakRegions_with_two_absent_chr_as_chrList <- function() {
    chrList <- Seqinfo(paste0("chr", c(32,1,40)), 
                        c(135534747, 249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                        Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                        peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                        chrInfo = chrList), error = conditionMessage)
    exp <- paste0("At least one chromosome name present in chrList is ",
                  "not present in peak : chr32, chr40")
    message <- paste0(" findConsensusPeakRegions_with_two_absent_chr_as",
                    "_chrList - Absent chromosomes in chrList did ", 
                    "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when zero as extendingSize
test.findConsensusPeakRegions_with_zero_as_extendingSize <- function() {
    chrList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], chrInfo = chrList,
            extendingSize = 0),
            error = conditionMessage)
    exp <- "extendingSize must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_with_zero_as_extendingSize",
                    " - Xero as extendingSize did not generated ", 
                    "expected error.")
    checkEquals(obs, exp, msg = message)
}
    
## Test the result when negative as extendingSize
test.findConsensusPeakRegions_with_negative_as_extendingSize <- function() {
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], chrInfo = testList, 
            extendingSize = -90),
            error = conditionMessage)
    exp <- "extendingSize must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_with_negative_as_extendingSize",
                    " - Negative as extendingSize did not generated ", 
                    "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when string as extendingSize
test.findConsensusPeakRegions_with_string_as_extendingSize <- function() {
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], chrInfo = testList,
            extendingSize = "444"), error = conditionMessage)
    exp <- "extendingSize must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_with_string_as_extendingSize",
                      " - String as extendingSize did not generated ", 
                      "expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when string as includeAllPeakRegion
test.findConsensusPeakRegions_string_as_includeAllPeakRegion <- function() {
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], chrInfo = testList,
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
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
            c(Hosa_A549_FOSL2_01_NarrowPeaks_partial, 
            Hosa_A549_FOXA1_01_NarrowPeaks_partial), 
            peaks = c(Hosa_A549_FOSL2_01_Peaks_partial, 
            Hosa_A549_FOXA1_01_Peaks_partial), chrInfo = testList,
            includeAllPeakRegion=333), error = conditionMessage)
    exp <- "includeAllPeakRegion must be a logical value"
    message <- paste0("findConsensusPeakRegions_numerical_as_",
                      "includeAllPeakRegion - Numerical as ", 
                      "includeAllPeakRegion did not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when string as minNbrExp
test.findConsensusPeakRegions_string_as_minNbrExp <- function() {
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                chrInfo = testList,
                minNbrExp = "444"), error = conditionMessage)
    exp <- "minNbrExp must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_string_as_minNbrExp",
                      " - String as minNbrExp did ", 
                      "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when zero as minNbrExp
test.findConsensusPeakRegions_zero_as_minNbrExp <- function() {
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2],  
                chrInfo = testList, 
                minNbrExp = 0), error = conditionMessage)
    exp <- "minNbrExp must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_zero_as_minNbrExp",
                      " - Zero as minNbrExp did ", 
                      "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when negative as minNbrExp
test.findConsensusPeakRegions_negative_as_minNbrExp <- function() {
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                                chrInfo = testList,
                                minNbrExp = -1), error = conditionMessage)
    exp <- "minNbrExp must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_negative_as_minNbrExp",
                    " - Negative as minNbrExp did ", 
                    "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when numerical as minNbrExp
test.findConsensusPeakRegions_numerical_as_minNbrExp <- function() {
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2],
                                chrInfo = testList,
                                minNbrExp = 9.3), error = conditionMessage)
    exp <- "minNbrExp must be a non-negative integer"
    message <- paste0(" findConsensusPeakRegions_numerical_as_minNbrExp",
                    " - Numerical as minNbrExp did ", 
                    "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when list of integers as minNbrExp
test.findConsensusPeakRegions_list_of_integers_as_minNbrExp <- function() {
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2],
                            chrInfo = testList,
                            minNbrExp= c(9L, 3L)), error = conditionMessage)
    exp <- "minNbrExp must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_list_of_integers_as_minNbrExp",
                    " - List of integers as minNbrExp did ", 
                    "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}


## Test the result when zero as nbrThreads
test.findConsensusPeakRegions_zero_as_nbrThreads<- function() {
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                                chrInfo = testList,
                                nbrThreads = 0), error = conditionMessage)
    exp <- "nbrThreads must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_zero_as_nbrThreads",
                    " - Zero as nbrThreads did ", 
                    "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when negative as nbrThreads
test.findConsensusPeakRegions_negative_as_nbrThreads <- function() {
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                                chrInfo = testList,
                                nbrThreads = -1), error = conditionMessage)
    exp <- "nbrThreads must be a non-negative integer"
    message <- paste0(" findConsensusPeakRegions_negative_as_nbrThreads",
                    " - Negative as nbrThreads did ", 
                    "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when numerical as nbrThreads
test.findConsensusPeakRegions_numerical_as_nbrThreads <- function() {
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                                Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                                peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                                chrInfo = testList, 
                                nbrThreads = 9.3), error = conditionMessage)
    exp <- "nbrThreads must be a non-negative integer"
    message <- paste0("findConsensusPeakRegions_numerical_as_nbrThreads",
                    " - Numerical as nbrThreads did ", 
                    "not generated expected error.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when list of integers as nbrThreads
test.findConsensusPeakRegions_list_of_integers_as_nbrThreads <- function() {
    testList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- tryCatch(findConsensusPeakRegions(narrowPeaks = 
                            Hosa_A549_FOSL2_01_NarrowPeaks_partial[1:2], 
                            peaks = Hosa_A549_FOSL2_01_Peaks_partial[1:2], 
                            chrInfo = testList, 
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

## Test the result when only one chromosome and one experiment
test.findConsensusPeakRegions_for_one_chromosome_one_experiment <- function() {
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
    testList <- Seqinfo(c("chr1"), c(249250621), NA)
    obs <- findConsensusPeakRegions(narrowPeaks = 
                c(Hosa_A549_FOSL2_01_NarrowPeaks_partial, 
                Hosa_A549_FOXA1_01_NarrowPeaks_partial), 
                peaks = c(Hosa_A549_FOSL2_01_Peaks_partial, 
                Hosa_A549_FOXA1_01_Peaks_partial), chrInfo = testList)
    message <- paste0(" findConsensusPeakRegions_for_one_chromosome_one",
                    "_experiment - When only one chromosome and one ", 
                    "experiment did not generated expected results.")
    checkEquals(obs$consensusRanges, exp, msg = message)
}

## Test the result when ALL as chrList
test.findConsensusPeakRegions_when_ALL <- function() {
    seqinfo <- Seqinfo(paste0("chr", c(1,10)), NA, NA, NA)
    exp <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(12,7)),
                   ranges = IRanges(start = c(249119914, 249120334, 249123074, 
                                              249132040, 249133011, 249152098, 
                                              249152823, 249153205, 249157198,
                                              249167214, 249167809, 249199968,
                                              179374,    182194,    183469,
                                              285046,    312979,    343055, 
                                              348698),
                                    end = c(249120424, 249121174, 249123574, 
                                            249132673, 249133517, 249152644, 
                                            249153397, 249153705, 249157698, 
                                            249167714, 249168802, 249200468,
                                            179874,    182694,    183969,
                                            285546,    313479,    343555, 
                                            349198)), seqinfo = seqinfo)
    chrList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- findConsensusPeakRegions(narrowPeaks = 
                    c(Hosa_A549_FOSL2_01_NarrowPeaks_partial, 
                    Hosa_A549_FOXA1_01_NarrowPeaks_partial), 
                    peaks = c(Hosa_A549_FOSL2_01_Peaks_partial, 
                    Hosa_A549_FOXA1_01_Peaks_partial), chrInfo = chrList)
    message <- paste0("findConsensusPeakRegions_for_one_chromosome ",
                      " - When \"ALL\" as chrList did", 
                      "not generated expected results.")
    checkEquals(obs$consensusRanges, exp, msg = message)
}

## Test the result when ALL as chrList and 2 as minNbrExp
test.findConsensusPeakRegions_when_ALL_with_minNbrExp_two <- function() {
    seqinfo <- Seqinfo(paste0("chr", c(1,10)), NA, NA, NA)
    exp <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(9, 3)),
                   ranges = IRanges(start = c(249119914, 249120334, 249123074, 
                                              249132040,  
                                              249152823, 249153205, 
                                              249167214, 249167809, 249199968,
                                              179374,    312979,    343055),
                                    end = c(249120424, 249121174, 249123574, 
                                            249132673,  
                                            249153397, 249153705,  
                                            249167714, 249168802, 249200468,
                                            179874,    313479,    343555)), 
                                    seqinfo = seqinfo)
    chrList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- findConsensusPeakRegions(narrowPeaks = 
                            c(Hosa_A549_FOSL2_01_NarrowPeaks_partial, 
                            Hosa_A549_FOXA1_01_NarrowPeaks_partial), 
                            peaks = c(Hosa_A549_FOSL2_01_Peaks_partial, 
                            Hosa_A549_FOXA1_01_Peaks_partial), 
                            chrInfo = chrList, minNbrExp = 2)
    message <- paste0("findConsensusPeakRegions_when_ALL_with_minNbrExp_two",
                      " - When \"ALL\" as chrList and two as minNbrExp did", 
                      "not generated expected results.")
    checkEquals(obs$consensusRanges, exp, msg = message)
}

## Test the result when ALL as chrList and 2 as minNbrExp and no expending region
test.findConsensusPeakRegions_ALL_with_minNbrExp_two_no_expending <- function() {
    seqinfo <- Seqinfo(paste0("chr", c(1,10)), NA, NA, NA)
    exp <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(9,3)),
                   ranges = IRanges(start = c(249119924, 249120334, 249123074, 
                                              249132040,  
                                              249152864, 249153205, 
                                              249167214, 249167809, 249199968,
                                              179374,    312979,    343055),
                                    end = c(249120424, 249120834, 249123574, 
                                            249132540,  
                                            249153364, 249153705,  
                                            249167714, 249168309, 249200468,
                                            179874,    313479,    343555)), 
                   seqinfo = seqinfo)
    chrList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- findConsensusPeakRegions(narrowPeaks = 
                            c(Hosa_A549_FOSL2_01_NarrowPeaks_partial, 
                            Hosa_A549_FOXA1_01_NarrowPeaks_partial), 
                            peaks = c(Hosa_A549_FOSL2_01_Peaks_partial, 
                            Hosa_A549_FOXA1_01_Peaks_partial), 
                            chrInfo = chrList,
                            minNbrExp = 2, includeAllPeakRegion = FALSE)
    message <- paste0("findConsensusPeakRegions_ALL_with_minNbrExp_two_no_",
                      "expending - When two as minNbrExp ",
                      "and no expending region did ", 
                      "not generated expected results.")
    checkEquals(obs$consensusRanges, exp, msg = message)
}

## Test the result when ALL as chrList and 2 as minNbrExp and no expending region
test.findConsensusPeakRegions_ALL_with_size_50_minNbrExp_two_no_expending <- function() {
    seqinfo <- Seqinfo(paste0("chr", c(1,10)), NA, NA, NA)
    exp <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(4,2)),
                   ranges = IRanges(start = c(249123274, 249167414, 249168009, 
                                            249200168,    179574,    343255),
                                    end = c(249123374, 249167514, 249168109, 
                                            249200268,    179674,    343355)), 
                   seqinfo = seqinfo)
    chrList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- findConsensusPeakRegions(narrowPeaks = 
                            c(Hosa_A549_FOSL2_01_NarrowPeaks_partial, 
                            Hosa_A549_FOXA1_01_NarrowPeaks_partial), 
                            peaks = c(Hosa_A549_FOSL2_01_Peaks_partial, 
                            Hosa_A549_FOXA1_01_Peaks_partial), 
                            chrInfo = chrList,
                            minNbrExp = 2, extendingSize = 50,
                            includeAllPeakRegion = FALSE)
    message <- paste0("findConsensusPeakRegions_ALL_with_size_50_minNbrExp_",
                      "two_no_expending - When \"ALL\" as chrList, two as ",
                      "minNbrExp, size of 50 and no expending region did ", 
                      "not generated expected results.")
    checkEquals(end(obs$consensusRanges)-start(obs$consensusRanges), 
                rep(100L, 6), msg = message)
    checkEquals(obs$consensusRanges, exp, msg = message)
}

## Test that left boundary inferior to zero is set to zero
test.findConsensusPeakRegions_ALL_with_one_as_left_boundary <- function() {
    seqinfo <- Seqinfo(paste0("chr", c(1,10)), NA, NA, NA)
    exp1Peak <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(1,1)),
                   ranges = IRanges(start = c(10, 40), end = c(10, 40)), 
                   name=c("peak1", "peak2"), 
                   seqinfo = seqinfo)
    names(exp1Peak)<-rep("exp1", 2)
    exp1NarrowPeak <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(1,1)),
                        ranges = IRanges(start = c(2, 34),
                                         end = c(33, 54)), 
                        name=c("peak1", "peak2"), 
                        seqinfo = seqinfo)
    names(exp1NarrowPeak)<-rep("exp1", 2)
    exp2Peak <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(1,1)),
                    ranges = IRanges(start = c(15, 35),
                                     end = c(15, 35)), 
                    name=c("peak1", "peak2"), 
                    seqinfo = seqinfo)
    names(exp2Peak)<-rep("exp2", 2)
    exp2NarrowPeak <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(1,1)),
                        ranges = IRanges(start = c(11, 32),
                                         end = c(19, 55)), 
                        name=c("peak1", "peak2"), 
                        seqinfo = seqinfo)
    names(exp2NarrowPeak)<-rep("exp2", 2)
    chrList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- findConsensusPeakRegions(narrowPeaks = 
                                        c(exp1NarrowPeak, 
                                          exp2NarrowPeak), 
                                    peaks = c(exp1Peak, 
                                              exp2Peak), 
                                    chrInfo = chrList,
                                    minNbrExp = 2, extendingSize = 100,
                                    includeAllPeakRegion = FALSE)
    exp <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(1,1)),
                   ranges = IRanges(start = c(1, 1),
                                    end = c(112, 137)), 
                   seqinfo = seqinfo)
    message <- paste0("findConsensusPeakRegions_ALL_with_one_as_left_bondary",
                      " - When left boubdary zero or negative, the boundary ", 
                      "is not modified to generate expected results.")
    checkEquals(obs$consensusRanges, exp, msg = message)
}

## Test that right boundary superior to chromosome length is set 
## to chromosome length
test.findConsensusPeakRegions_ALL_with_superior_right_boundary <- function() {
    seqinfo <- Seqinfo(paste0("chr", c(1,10)), NA, NA, NA)
    exp1Peak <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(1,1)),
                        ranges = IRanges(start = c(249250617, 135534737), 
                                        end = c(249250617, 135534737)), 
                        name=c("peak1", "peak2"), 
                        seqinfo = seqinfo)
    names(exp1Peak)<-rep("exp1", 2)
    exp1NarrowPeak <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(1,1)),
                              ranges = IRanges(start = c(249250614, 135534717),
                                               end = c(249250619, 135534737)), 
                              name=c("peak1", "peak2"), 
                              seqinfo = seqinfo)
    names(exp1NarrowPeak)<-rep("exp1", 2)
    exp2Peak <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(1,1)),
                        ranges = IRanges(start = c(249250619, 135534740),
                                         end = c(249250619, 135534740)), 
                        name=c("peak1", "peak2"), 
                        seqinfo = seqinfo)
    names(exp2Peak)<-rep("exp2", 2)
    exp2NarrowPeak <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(1,1)),
                              ranges = IRanges(start = c(249250600, 135534710),
                                               end = c(249250620, 135534746)), 
                              name=c("peak1", "peak2"), 
                              seqinfo = seqinfo)
    names(exp2NarrowPeak)<-rep("exp2", 2)
    chrList <- Seqinfo(paste0("chr", c(1,10)), c(249250621, 135534747), NA)
    obs <- findConsensusPeakRegions(narrowPeaks = 
                                        c(exp1NarrowPeak, 
                                          exp2NarrowPeak), 
                                    peaks = c(exp1Peak, 
                                                exp2Peak), 
                                    chrInfo = chrList,
                                    minNbrExp = 2, extendingSize = 100,
                                    includeAllPeakRegion = FALSE)
    exp <- GRanges(seqnames = Rle(c("chr1", "chr10"),c(1,1)),
                   ranges = IRanges(start = c(249250518, 135534638),
                                    end = c(249250621, 135534747)), 
                   seqinfo = seqinfo)
    message <- paste0("findConsensusPeakRegions_ALL_with_one_as_left_bondary",
                      " - When right boubdary superior to chromosome length ", 
                      "the boundary ", 
                      "is not modified to generate expected results.")
    checkEquals(obs$consensusRanges, exp, msg = message)
}