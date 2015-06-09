###################################################
# Created by Astrid Louise Deschenes
# 2015-06-01
###################################################

###################################################
## Test the findConsensusPeakRegionsIntern.R functions
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "consensusSeekeR" )
}

### }}}

###################################################
## Test the isInteger function
###################################################

## Test the result when a integer is passed to the function
test.isInteger_with_integer <- function() {
    obs <- tryCatch(consensusSeekeR:::isInteger(value = 444), 
                    error=conditionMessage)
    message <- paste0("isInteger_with_integer() ",
                    "- A integer did not generated TRUE.")
    checkEquals(obs, TRUE, msg = message)
}

## Test the result when a string is passed to the function
test.isInteger_with_string <- function() {
    obs <- tryCatch(consensusSeekeR:::isInteger(value = "444"), 
                    error=conditionMessage)
    message <- paste0("isInteger_with_integer() ",
                    "- A string did not generated FALSE.")
    checkEquals(obs, FALSE, msg = message)
}

## Test the result when a list is passed to the function
test.isInteger_with_list_of_integers <- function() {
    obs <- tryCatch(consensusSeekeR:::isInteger(value = list(a=444, b=333)), 
                    error=conditionMessage)
    message <- paste0("isInteger_with_integer() ",
                    "- A list did not generated FALSE.")
    checkEquals(obs, FALSE, msg = message)
}

## Test the result when a vector is passed to the function
test.isInteger_with_vector_of_integers <- function() {
    obs <- tryCatch(consensusSeekeR:::isInteger(value = c(444,333)), 
                    error=conditionMessage)
    message <- paste0("isInteger_with_integer() ",
                "- A vector of integers did not generated FALSE.")
    checkEquals(obs, FALSE, msg = message)
}

