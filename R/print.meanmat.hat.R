print.meanmat.hat <-
function (x, ...) 
{
 cat("ESTIMATION OF THE MEAN MATRIX", "\n")
 cat("Sample Size           = ",x$N,"\n")
 cat("Row Variables         = ",x$n.rows,"\n")
 cat("Columns Variables     = ",x$n.cols,"\n")
 cat("\nEstimated Mean Matrix [1:5,1:5] =\n"); print(round(x$estmeanmat[1:min(5,x$n.rows),1:min(5,x$n.cols)],4))
}
