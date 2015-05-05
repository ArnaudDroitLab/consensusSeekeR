
#######################################################
## Load libraries
#######################################################
# 
# library(rtracklayer)
# library(GenomicRanges)
# library(XVector)
# library(stringr)
# library(BiocParallel)
# library(BiocGenerics)

#######################################################
## Load data
#######################################################

# narrowPeakSMC1 <- import("..//..//data//data_case_02//HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1_R1.trim.hg19.sorted.bam_peaks.narrowPeak.net.bed",format="bed",genome="hg19")
# narrowPeakMED1 <- import("..//..//data//data_case_02//HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1_R1.trim.hg19.sorted.bam_peaks.narrowPeak.net.bed",format="bed",genome="hg19")
# narrowPeakNIPBL <- import("..//..//data//data_case_02//HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1_R1.trim.hg19.sorted.bam_peaks.narrowPeak.net.bed",format="bed",genome="hg19")
# peakSMC1 <- import("..//..//data//data_case_02//HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1_R1.trim.hg19.sorted.bam_summits.net.bed",format="bed",genome="hg19")
# peakMED1 <- import("..//..//data//data_case_02//HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1_R1.trim.hg19.sorted.bam_summits.net.bed",format="bed",genome="hg19")
# peakNIPBL <- import("..//..//data//data_case_02//HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1_R1.trim.hg19.sorted.bam_summits.net.bed",format="bed",genome="hg19")


#######################################################
## Group peak data together
#######################################################
# 
# allPeaks <- sort(c(peakSMC1, peakMED1, peakNIPBL))
# allNarrowPeaks <- sort(c(narrowPeakMED1, narrowPeakNIPBL, narrowPeakSMC1))
# 
# rm(narrowPeakMED1, narrowPeakNIPBL, narrowPeakSMC1, peakMED1, peakNIPBL, peakSMC1)

#######################################################
## Load data
#######################################################
findCommonRegions <- function(chrName, allPeaks, allNarrowPeaks) {
    ## Keep data related to the selected chromosome
    peaks <- subset(allPeaks, seqnames(allPeaks) == chrName)
    narrowPeaks <- subset(allNarrowPeaks, seqnames(allNarrowPeaks) == chrName)
    
    regions <- GRanges()
    strangeRegions <- GRanges()
    current <- NULL
    rightBoundary <- NULL
    namesVec <- vector()
    bad <- FALSE
    pos <- 1
    
    repeat  {
        current <- peaks[pos]
        rightBoundaryNew <- start(current)
        rightBoundary <- NULL
        set <- NULL
        setNew <- NULL
        bad <- FALSE
        repeat {
            set <- setNew
            rightBoundary <- rightBoundaryNew
            # Find peaks that overlaps the selected 1000 bases region
            overlaps <- findOverlaps(query = GRanges(seqnames = c(as.character(seqnames(current))),
                                                     ranges=c(IRanges(rightBoundary, rightBoundary+1000))),
                                     subject = peaks)
            setNew <- peaks[subjectHits(overlaps)]
            if (!(current$name %in% setNew$name)) {
                # The current peak is not included in the current region
                # The region will not be selected
                strangeRegions <- append(strangeRegions, current)
                bad <- TRUE
                break
            }
            # Use the median of the peaks to set the new right boundary
            rightBoundaryNew <- median(start(setNew))-500
            # Stop loop when the overlaping peaks are stable or when no peaks are found
            if (!is.null(set) && (length(set) == length(setNew)) && all(set == setNew)) break
        }
        
        if (bad) {
            # Treat the next position
            pos<-pos + 1
        } else {
            # Keep region only when peaks from more than one experience are present
            browser()
            short_names <- sapply(X = set$name, function(x) stringr::str_split(string = x, pattern = ".bam")[[1]][1])
            if (length(unique(short_names)) > 1) {
                # Create one final region using the narrow information for each peak present
                minPos <- ifelse(min(BiocGenerics::start(narrowPeaks[narrowPeaks$name %in%  set$name])), rightBoundary)
                maxPos <- rightBoundary + 1000
                for (i in unique(short_names)) {
                    peaksForOneExp <- set[short_names == i]
                    
                    firstPeak <- peaksForOneExp[1]
                    lastPeak <- peaksForOneExp[length(peaksForOneExp)]
                    
                    newMax <- BiocGenerics::end(narrowPeaks[narrowPeaks$name %in%  firstPeak$name])
                    maxPos <- ifelse(newMax > maxPos, newMax, maxPos)
                    
                    newMin <- BiocGenerics::start(narrowPeaks[narrowPeaks$name %in%  lastPeak$name])
                    minPos <- ifelse(newMin < minPos, newMin, minPos)
                }
                
                newRegion <- GRanges(seqnames = as.character(seqnames(set[1])), IRanges(minPos, maxPos))
                regions <- append(regions, newRegion)
                # Update overlapping peaks
                overlaps <- findOverlaps(query = newRegion, subject = peaks)
            }
            # Treat the position following last peak present in new region
            
            pos<-max(subjectHits(overlaps))+1
        }
        # Stop loop when all peaks are treated
        if (pos >= length(peaks)) break
        #if (pos >= 10) break
    }
    return(list(regions, strangeRegions))
}

# results<-bplapply(levels(seqnames(allPeaks))[1:24], FUN = findCommonRegions,
#                   allPeaks=allPeaks, allNarrowPeaks=allNarrowPeaks,
#                   BPPARAM = MulticoreParam(workers = 13))
# 
# finalRegions<-GRanges()
# strangeRegions<-GRanges()
# for (i in 1:length(results)) {
#     finalRegions<-c(finalRegions, results[[i]][[1]])
#     strangeRegions<-c(strangeRegions, results[[i]][[2]])
# }
# 
# finalRegions<-sort(finalRegions)
# strangeRegions <- sort(strangeRegions)
# 
# save(finalRegions, file="finalRegions_CASE_02_Method_02_NEW.Rda")
# save(strangeRegions, file="strangeRegions_CASE_02_Method_02_NEW.Rda")
# export(finalRegions, file("finalRegions_CASE_02_Method_02_NEW.bed"))
# export(strangeRegions, file("strangeRegions_CASE_02_Method_02_NEW.bed"))
