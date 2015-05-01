#' @title TODO
#' 
#' @description TODO
#' 
#' @param peaksBEDlist an object of \code{class} "formula" which contains a symbolic
#'          model formula.
#' @param narrowpeaksBEDlist a \code{data.frame} containing the variables in the model.
#' @param chrList a \code{GRanges} indicating which column from the 
#'          \code{data} must be added to the formula. When \code{NULL}, no
#'          new term is added. Default : \code{NULL}.
#' @param padding a \code{GRanges}. Default = 500.
#' @param minNbrExp a \code{} Default = 1.
#' 
#' @return an object of \code{class} "commonFeatures". 
#' 
#' @author Astrid Louise Deschenes
#' @importFrom BiocGenerics start end
#' @importFrom stringr str_split
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @keywords internal
findCommonFeatures <- function(peaksBEDlist, narrowpeaksBEDlist, chrList, 
                               padding = 500, minNbrExp = 1) {