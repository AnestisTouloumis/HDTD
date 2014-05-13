print.meanmat.ts <-
function (x, ...) 
{
 cat("MEAN MATRIX TEST", "\n")
 cat("Sample Size       =",x$N,"\n")
 cat("Row Variables     =",x$n.rows,"\n")
 cat("Columns Variables =",x$n.cols,"\n")
 cat("\nHypothesis Test\n")
 cat("H_0:",x$n.groups,"prespecified group(s) of",x$voi,"with the same mean vector\n")
 cat("vs.\n")
 cat("H_1: not H_0 \n")
 if(x$n.groups==2) cat("\nThe number of",x$voi,"in each of the two predefined groups are", paste(x$group.sizes,collapse=" and "), ".\n")
 if(x$n.groups>2) cat("\nThe number of",x$voi,"in the",x$n.groups,"predefined groups are",paste(c(paste(x$group.sizes[1:(x$n.groups-1)],collapse=", "),x$group.sizes[x$n.groups]),collapse=" and ")," respectively.\n")
 cat("\nTest Statistic =",round(x$statistic,4),", p-value =",ifelse(x$p.value<=0.0001,"<0.0001",round(x$p.value,4)),"\n")
}
