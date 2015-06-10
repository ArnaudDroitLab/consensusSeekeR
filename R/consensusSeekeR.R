#' consensusSeekeR: Detection of consensus peak regions inside a group of
#' experiments using narrowPeak files
#'
#' This package compares multiple narrowPeak files to extract common
#' consensus peak regions. The size of the analyzed region is adjustable
#' as well as the number of narrowPeak files in which a peak must be detected
#' in the potential region to mark the region as a consensus peak region.
#'
#' @docType package
#' @name consensusSeekeR-package
#' @aliases consensusSeekeR-package consensusSeekeR
#' @author  Astrid Louise Deschenes,
#'  Fabien Claude Lamaze,
#'  Pascal Belleau and
#'  Arnaud Droit
#'
#' Maintainer: Astrid Louise Deschenes <astrid-louise.deschenes@@crchudequebec.ulaval.ca>
#' @keywords package
NULL

#' Regions of sites with the greatest evidence of transcription factor binding
#' for the FOSL2 transcription factor (for demonstration purpose)
#'
#' Regions of sites with the greatest evidence of transcription factor binding
#' for the
#' FOSL2 transcription factor (DCC accession: ENCFF000MZT), calculated
#' using the MACS peak caller (Zhang et al., 2008),
#' for regions chr1:249120200-249250621 and chr10:1-370100
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#' This dataset is associated to the
#'
#' @name Hosa_A549_FOSL2_01_NarrowPeaks
#' @docType data
#' @format A \code{GRanges} containing a metadata field called "name"
#'      which identified the peak associated to the narrow peak region.
#'      Each \code{GRanges} entry has a
#'      different peak name.
#' @usage data(Hosa_A549_FOSL2_01_NarrowPeaks)
#' @seealso
#'  \itemize{
#'    \item \code{\link{Hosa_A549_FOSL2_01_Peaks}} { the associate peaks
#'                  dataset}
#'    \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#'}
#' @keywords datasets
NULL

#' Sites with the greatest evidence of transcription factor binding for the
#' FOSL2 transcription factor (for demonstration purpose)
#'
#' Sites with the greatest evidence of transcription factor binding for the
#' FOSL2 transcription factor (DCC accession: ENCFF000MZT), calculated
#' using the MACS peak caller (Zhang et al., 2008),
#' for regions chr1:249120200-249250621 and chr10:1-370100
#' from
#' the Encyclopedia of DNA Elements (ENCODE) data (Dunham I et al. 2012).
#'
#' @name Hosa_A549_FOSL2_01_Peaks
#' @docType data
#' @format A \code{GRanges} containing a metadata field called "name"
#'      which identified the peak. Each \code{GRanges} entry has a
#'      different peak name.
#' @usage data(Hosa_A549_FOSL2_01_Peaks)
#' @seealso
#'  \itemize{
#'    \item \code{\link{Hosa_A549_FOSL2_01_NarrowPeaks}} { the associate
#'                  narrow peaks dataset}
#'    \item \code{\link{findConsensusPeakRegions}} {for extracting regions
#'                  sharing the same features in more than one experiment.}
#' }
#' @references
#'  \itemize{
#'  \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia of
#'  DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#'  \item Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, et
#'  al. Model-based analysis of ChIP-Seq (MACS). Genome  Biol. 2008, 9(9):R137.
#'  }
#' @keywords datasets
NULL

#' TODO Data Set
#'
#' TODO
#'
#' @name Hosa_A549_FOSL2_01_NarrowPeaks_partial
#' @docType data
#' @format A \code{GRanges}
#' @references
#'  \itemize{
#'  \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia of DNA
#'  elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#'  \item Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, et al.
#'  Model-based analysis of ChIP-Seq (MACS). Genome  Biol. 2008, 9(9):R137.
#'  }
#' @keywords datasets
NULL

#' TODO Data Set
#'
#' TODO
#'
#' @name Hosa_A549_FOSL2_01_Peaks_partial
#' @docType data
#' @format A \code{GRanges}
#' @usage data(Hosa_A549_FOSL2_01_Peaks_partial)
#' @references
#'  \itemize{
#'  \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia of DNA
#'  elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#'  \item Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, et al.
#'  Model-based analysis of ChIP-Seq (MACS). Genome  Biol. 2008, 9(9):R137.
#'  }
#' @keywords datasets
NULL

#' TODO Data Set
#'
#' TODO
#'
#' @name Hosa_A549_FOXA1_01_NarrowPeaks_partial
#' @docType data
#' @format A \code{GRanges}
#' @usage data(Hosa_A549_FOXA1_01_NarrowPeaks_partial)
#' @references
#'  \itemize{
#'  \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia of DNA
#'  elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#'  \item Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, et al.
#'  Model-based analysis of ChIP-Seq (MACS). Genome  Biol. 2008, 9(9):R137.
#'  }
#' @keywords datasets
NULL

#' TODO Data Set
#'
#' @name Hosa_A549_FOXA1_01_Peaks_partial
#' @docType data
#' @format A \code{GRanges}
#' @usage data(Hosa_A549_FOXA1_01_Peaks_partial)
#' @references
#'  \itemize{
#'  \item Dunham I, Kundaje A, Aldred SF, et al. An integrated encyclopedia of DNA
#'  elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74.
#'  \item Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, et al.
#'  Model-based analysis of ChIP-Seq (MACS). Genome  Biol. 2008, 9(9):R137.
#'  }
#' @keywords datasets
NULL