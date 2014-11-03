\name{transposedata}
\alias{transposedata}
\title{
Interchanging the Row and Column Variables in Transposable Data
}
\description{
This function interchanges the row and column variables in transposable data so that the original row variables will be treated as column variables and the original column variables as row variables.
}
\usage{
transposedata(datamat, N)
}
\arguments{
  \item{datamat}{numeric matrix containing the transposable data.
}
  \item{N}{positive integer number indicating the sample size, i.e., the number of subjects.
}
}
\details{
It is assumed that there are \code{nrow(datamat)} row variables and \code{ncol(datamat)}/\code{N} column variables in \code{datamat}. Further, \code{datamat} should be written in such a way that every \code{ncol(datamat)}/\code{N} consecutive columns belong to the same subject and the order of the column variables in each block is preserved across subjects.
}
\value{
Returns a matrix with \code{ncol(datamat)} rows and \code{nrow(datamat)}\code{N} columns. 
}
\author{
Anestis Touloumis
}
\seealso{
\code{\link{centerdata}} and \code{\link{orderdata}}.
}
\examples{
data(VEGFmouse)
## Transposing the VEGF dataset.
VEGFtr <- transposedata(VEGFmouse,40)
}

