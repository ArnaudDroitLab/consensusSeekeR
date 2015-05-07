#' @title Extract narrow regions and peaks from narrrowPeak file
#' 
#' @description Read a narrowPeak file and extract the narrow regions and/or
#'          the peaks, as specified by used. The narrowPeak file must fit the 
#'          UCSC specifications. See 
#'          \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format12} for more
#'          details.
#' 
#' @param file_path the name of the file.
#' @param extractRegions a \code{logical} indicating if the narrow regions must
#'          be extracted. If \code{TRUE}, a \code{GRanges} containing the narrow
#'          regions will be returned. Default = \code{TRUE}.
#' @param extractPeaks a \code{logical} indicating if the peaks must
#'          be extracted. If \code{TRUE}, a \code{GRanges} containing the peaks
#'          will be returned. Default = \code{TRUE}.
#' 
#' @return a \code{list} containing 2 entries:
#'      \itemize{
#'          \item narrowPeak a {\code{GRanges}} containing 
#'              the narrow regions extracted from the file.
#'          \item peak a {\code{GRanges}} containing 
#'              the peaks extracted from the file.
#'      }
#' 
#' @author Astrid Louise Deschenes
#' @importFrom BiocGenerics start end
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @export
readNarrowPeak<- function(file_path, extractRegions = TRUE, 
                                extractPeaks = TRUE) {
    # Parameters validation
    if (!file.exists(file_path)) {
        stop("No such file: ", file_path)
    }
    
    if (!is.logical(extractRegions)) {
        stop("extractRegions must be a logical value")
    }
    
    if (!is.logical(extractPeaks)) {
        stop("extractPeaks must be a logical value")
    }
    
    if (!extractRegions && !extractPeaks) {
        stop("extractPeaks and extractRegions cannot be both FALSE")
    }
    
    chunk_size = 250
    chunk = scan(file_path, what="character", nmax=chunk_size, 
                 strip.white=TRUE, sep="\n", quiet=TRUE)
    skip_n = suppressWarnings(min(grep("^chr(\\d+|\\w+)\\s+\\d+",chunk)) - 1)
    
    if (is.infinite(skip_n)) {
        stop("No valid chromosomes detected within first 250 ", 
                "lines of BED file \"", file_path , "\"")
    }
    
    peaks = read.table(file_path, header=FALSE, skip=skip_n)
    peaks = peaks[,1:10];
    names(peaks) = c("chrom","start", "end", "name", "score", "strand", 
                        "signalValue", "pValue", "qValue", "peak")
    
    if (any(peaks$start < 0) || any(peaks$end < 0)) {
        stop("Start ans end positions of peaks should all be >= 0.")
    }
    
    regionResult = list();
    peakResult = list();
#     for (chr in unique(peaks$chrom)) {
#         peaks_chrom = subset(peaks, peakschrom == chr)
#         chromResult[[chr]] <- GRanges(seqnames = as.character(chr), 
#                             IRanges(start=(peaks_chrom$start + 1L), 
#                             end=peaks_chrom$end),
#                             name = as.character(peaks_chrom$name))
#         peakResult[[chr]] <- GRanges(seqnames = as.character(chr), 
#                             IRanges(start=(peaks_chrom$start + 1L + 
#                             peaks_chrom$peak), end=(peaks_chrom$start + 1L + 
#                             peaks_chrom$peak)), 
#                             name = as.character(peaks_chrom$name))
#     }

    if (extractRegions) {
        regionResult <- GRanges(seqnames = as.character(peaks$chrom), 
                            IRanges(start=(peaks$start + 1L), 
                            end=peaks$end),
                            name = as.character(peaks$name))
    }

    if (extractPeaks) {
        peakResult <- GRanges(seqnames = as.character(peaks$chrom), 
                            IRanges(start=(peaks$start + 1L + 
                            peaks$peak), end=(peaks$start + 1L + 
                            peaks$peak)), 
                            name = as.character(peaks$name))
    }

    return(list(narrowPeak = regionResult, peak = peakResult))
}