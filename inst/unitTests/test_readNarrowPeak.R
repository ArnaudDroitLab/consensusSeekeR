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
