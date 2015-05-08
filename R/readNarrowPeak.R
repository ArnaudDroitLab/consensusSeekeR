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
        stop("No such file \"", file_path, "\"")
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
    
    data_size = 250
    data = scan(file_path, what="character", nmax=chunk_size, 
                 strip.white=TRUE, sep="\n", quiet=TRUE)
    
    grepRes <- grep(paste0("^\\S+(\\s+\\d+){2}\\s+\\S+\\s+\\d+\\s+[-\\+*\\.]",
                "(\\s+[0-9\\.]+){3}\\s+\\d+"), data)
    
    if (length(grepRes) == 0) {
        stop("No valid chromosome detected within first ", data_size, 
             " lines of BED file \"", file_path , "\"")
    }
    skip_lines = min(grepRes) - 1
    
    peaks = read.table(file_path, header = FALSE, skip = skip_lines)
    peaks = peaks[,1:10];
    names(peaks) = c("chrom","start", "end", "name", "score", "strand", 
                        "signalValue", "pValue", "qValue", "peak")
    
    if (any(peaks$start < 0) || any(peaks$end < 0)) {
        stop("start and end positions of peaks should all be >= 0.")
    }
    
    if (any(levels(peaks$strand) == ".")) {
        levels(peaks$strand)[levels(peaks$strand) == "."] <- "*"
    }
    
    regionResult = list();
    peakResult = list();

    if (extractRegions) {
        regionResult <- GRanges(seqnames = as.character(peaks$chrom), 
                            IRanges(start=(peaks$start + 1L), 
                            end=peaks$end),
                            name = as.character(peaks$name),
                            score = as.integer(peaks$score),
                            signalValue = as.numeric(peaks$signalValue),
                            strand = Rle(as.character(peaks$strand)),
                            pValue = as.numeric(peaks$pValue),
                            qValue = as.numeric(peaks$qValue),
                            peak = as.integer(peaks$peak)
                            )
    }

    if (extractPeaks) {
        peakResult <- GRanges(seqnames = as.character(peaks$chrom), 
                            IRanges(start=(peaks$start + 1L + 
                            peaks$peak), end=(peaks$start + 1L + 
                            peaks$peak)), 
                            name = as.character(peaks$name),
                            score = as.integer(peaks$score),
                            signalValue = as.numeric(peaks$signalValue),
                            strand = Rle(as.character(peaks$strand)),
                            pValue = as.numeric(peaks$pValue),
                            qValue = as.numeric(peaks$qValue),
                            peak = as.integer(peaks$peak)
                            )
    }

    return(list(narrowPeak = regionResult, peak = peakResult))
}