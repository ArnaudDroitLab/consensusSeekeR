#' @title TODO
#' 
#' @description TODO
#' 
#' @param file_path a \code{vector} 
#' 
#' @return \code{0}
#' 
#' @author Astrid Louise Deschenes
#' @importFrom BiocGenerics start end
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @keywords internal
readNarrowPeak<- function(file_path) {
    if (!file.exists(file_path)) {
        stop("No such file: ", file_path)
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
    
    chromResult = list();
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
    chromResult<- GRanges(seqnames = as.character(peaks$chrom), 
                            IRanges(start=(peaks$start + 1L), 
                            end=peaks$end),
                            name = as.character(peaks$name))
    peakResult <- GRanges(seqnames = as.character(peaks$chrom), 
                            IRanges(start=(peaks$start + 1L + 
                            peaks$peak), end=(peaks$start + 1L + 
                            peaks$peak)), 
                            name = as.character(peaks$name))

    return(list(narrowPeak=chromResult, peak=peakResult))
}