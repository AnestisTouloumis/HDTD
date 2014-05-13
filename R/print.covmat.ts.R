print.covmat.ts <-
function (x, ...) 
{
 cat("SPHERICITY AND IDENTITY TESTS FOR THE ROW OR COLUMN VARIABLES", "\n")
 cat("Sample Size           = ",x$N,"\n")
 cat("Row Variables         = ",x$n.rows,"\n")
 cat("Columns Variables     = ",x$n.cols,"\n")
 cat("Variables Tested      = ",x$variables,"\n")
 cat("Centered Data         = ",x$centered,"\n")
 cat("\nSphericity test for the covariance matrix of the",x$variables,"\n")
 cat("Test Statistic =",round(x$sphericity.ts$statistic,4),", p-value =",ifelse(x$sphericity.ts$p.value<=0.0001,"<0.0001",round(x$sphericity.ts$p.value,4)),"\n")
 cat("\nIdentity Test for the covariance matrix of the",x$variables,"\n")
 cat("Test Statistic =",round(x$identity.ts$statistic,4),", p-value =",ifelse(x$identity.ts$p.value<=0.0001,"<0.0001",round(x$identity.ts$p.value,4)),"\n")
}
