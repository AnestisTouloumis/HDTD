print.covmat.hat <-
function (x, ...) 
{
 cat("ESTIMATION OF THE ROW AND/OR THE COLUMN COVARIANCE MATRIX", "\n")
 cat("Sample Size           = ",x$N,"\n")
 cat("Row Variables         = ",x$n.rows,"\n")
 cat("Column Variables      = ",x$n.cols,"\n")
 cat("Shrinking             = ",x$shrink,"\n")
 cat("Centered Data         = ",x$centered,"\n")
 cat("\nROW VARIABLES\n")
 if(!is.null(x$rows.intensity)) cat("Estimated Shrinkage Intensity =",round(x$rows.intensity,4),"\n")
 cat("Estimated Covariance Matrix [1:5,1:5] =\n"); print(round(x$rows.covmat[1:min(5,x$n.rows),1:min(5,x$n.rows)],4))
 cat("\nCOLUMN VARIABLES\n")
 if(!is.null(x$cols.intensity)) cat("Estimated Shrinkage Intensity =",round(x$cols.intensity,4),"\n")
 cat("Estimated Covariance Matrix [1:5,1:5] =\n"); print(round(x$cols.covmat[1:min(5,x$n.cols),1:min(5,x$n.cols)],4))
}
