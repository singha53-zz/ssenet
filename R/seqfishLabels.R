#' seqfishLabels
#'
#' Spatial cluster labels and SVM learned cell types for seqFISH
#'
#' @docType data
#'
#' @usage data(seqfishLabels)
#'
#' @format An object of class \code{"data.frame"};
#' \describe{
#' \item{seqfishLabels}{1597 cells x 5 variables (Cell ID, Spatial cluster, cell-type class, ID of cell type class and probability of the predicted cell-type label)}
#' }
#'
#' @keywords scRNAseq, seqFISH,
#'
#' @references Zhu \emph{et al.} Identification of spatially associated subpopulations by combining scRNAseq and sequential fluorescence in situ hybridization data. \emph{Nat Biotechnol.} 2018 Oct 29.
#' (\href{https://www.nature.com/articles/nbt.4260}{Zhu \emph{et al} 2018})
#'
#'
#' @examples
#' library(ssenet)
#' data(seqfishLabels)
"seqfishLabels"
