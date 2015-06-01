###################################################
# Created by Astrid Louise Deschenes
# 2015-05-11
###################################################

###################################################
## Test the readnarrowPeak() function
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "sharedBed" )
}

### }}}

file_FOSL2_Rep01 <- dir(system.file("extdata", package = "sharedBed"), 
            pattern = "Hosa_A549_FOSL2_ENCSR000BQO_ENCFF000MZT_peaks_part_chr1_and_chr12.narrowPeak$",
            full.names=TRUE)


## Test the result when the file doesn't exist
test.readNarrowPeak_no_existing_file <- function() {
    file <- "TiiTooTiiToo.narrowPeak"
    obs <- tryCatch(readNarrowPeak(file_path = file), 
                    error=conditionMessage)
    exp <- paste0("No such file \"", file, "\"")
    message <- paste0("readNarrowPeak_no_existing_file() ",
        "- A not existing file did not generated the expected exception.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when extractPeaks and extractRegions parameters
## are both FALSE
test.readNarrowPeak_both_parameters_FALSE <- function() {
    obs <- tryCatch(readNarrowPeak(file_path = file_FOSL2_Rep01, 
                    extractRegions = FALSE, extractPeaks = FALSE), 
                    error=conditionMessage)
    exp <- "extractPeaks and extractRegions cannot be both FALSE"
    message <- paste0("readNarrowPeak_no_existing_file() ",
                     "- Both extractPeaks and extractRegions set to ",
                    "FALSE did not generated the expected exception.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when extractPeaks is a string
test.readNarrowPeak_a_string_as_extractPeaks <- function() {
    obs <- tryCatch(readNarrowPeak(file_path = file_FOSL2_Rep01, 
                                    extractRegions = FALSE, 
                                    extractPeaks = "extremeYoga"), 
                                    error=conditionMessage)
    exp <- "extractPeaks must be a logical value"
    message <- paste0("readNarrowPeak_no_existing_file() ",
                      "- A string used as extractPeaks parameter ",
                      "did not generated the expected exception.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when extractPeaks is a numerical
test.readNarrowPeak_a_numerical_as_extractPeaks <- function() {
    obs <- tryCatch(readNarrowPeak(file_path = file_FOSL2_Rep01, 
                                   extractRegions = FALSE, 
                                   extractPeaks = 0.001), 
                    error=conditionMessage)
    exp <- "extractPeaks must be a logical value"
    message <- paste0("readNarrowPeak_no_existing_file() ",
                      "- A string used as extractPeaks parameter ",
                      "did not generated the expected exception.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when extractRegions is a string
test.readNarrowPeak_a_string_as_extractRegions <- function() {
    obs <- tryCatch(readNarrowPeak(file_path = file_FOSL2_Rep01, 
                                   extractRegions = "extremeReading", 
                                   extractPeaks = TRUE), 
                    error=conditionMessage)
    exp <- "extractRegions must be a logical value"
    message <- paste0("readNarrowPeak_no_existing_file() ",
                      "- A string used as extractRegions parameter ",
                      "did not generated the expected exception.")
    checkEquals(obs, exp, msg = message)
}

# Test the result when extractRegions is an numerical 
test.readNarrowPeak_a_numerical_as_extractRegions <- function() {
    obs <- tryCatch(readNarrowPeak(file_path = file_FOSL2_Rep01, 
                                   extractRegions = 0.11, 
                                   extractPeaks = TRUE), 
                    error=conditionMessage)
    exp <- "extractRegions must be a logical value"
    message <- paste0("readNarrowPeak_no_existing_file() ",
                      "- A numerical used as extractRegions parameter ",
                      "did not generated the expected exception.")
    checkEquals(obs, exp, msg = message)
}
