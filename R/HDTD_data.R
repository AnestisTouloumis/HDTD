#' Vascular Endothelial Growth Factor Mouse Dataset
#' 
#' Log2 normalized mouse gene expression data in the vascular endothelial
#' growth factor signalling pathway across multiple tissues.
#' 
#' 
#' @name VEGFmouse
#' @docType data
#' @format A data frame with 46 rows and 360 columns. The rows corresponds to
#' 46 genes in the VEGF signalling pathway. The column names indicate the mouse
#' and the tissue on which gene expression levels were measured.  Since there
#' are 40 mice and 9 tissues, we have a total of 360 columns. Every 9
#' consecutive columns belong to the same mouse and the tissues are ordered in
#' the same way in each mouse.
#' @source Zahn et al. (2007). AGEMAP: A gene expression database for aging in
#' mice. \emph{PLoS Genetics} \bold{3}, e201.
#' @keywords datasets
#' @examples
#' data(VEGFmouse)
#' ## Check the order of the tissues from the first mouse.
#' colnames(VEGFmouse[,1:9]) 
NULL
