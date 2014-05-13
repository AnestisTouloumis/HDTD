\name{covmat.ts}
\alias{covmat.ts}
\alias{print.covmat.ts}
\title{
Nonparametric Tests for the Row or Column Covariance Matrix
}
\description{
Testing the sphericity and identity hypotheses for the row or column covariance matrix.
}
\usage{
covmat.ts(datamat, N, voi = "rows", centered = FALSE)
}
\arguments{
  \item{datamat}{numeric matrix containing the transposable data.
}
 \item{N}{
positive integer number indicating the sample size, i.e., the number of subjects.
}
  \item{voi}{
character indicating if the test should be applied on the row or column covariance matrix. Options include "\code{rows}" or "\code{columns}".
}
  \item{centered}{
logical indicating if the transposable data are centered. Options include \code{TRUE} or \code{FALSE}.
}
}
\details{
It is assumed that there are \code{nrow(datamat)} row variables and \code{ncol(datamat)}/\code{N} column variables in \code{datamat}. Further, \code{datamat} should be written in such a way that every \code{ncol(datamat)}/\code{N} consecutive columns belong to the same subject and the order of the column variables in each block is preserved across subjects.

The tests are nonparametric and thus robust to departures from the matrix-variate normal model. 
}
\value{
It returns a list with components:
\item{sphericity.ts}{a list containing the test statistic and p-value of the sphericity test.}
\item{identity.ts}{a list containing the test statistic and p-value of the identity test.}
\item{N}{the sample size.}
\item{n.rows}{the number of row variables.}
\item{n.cols}{the number of column variables.}
\item{variables}{character indicating if the tests were applied to the row or column covariance matrix.}
\item{centered}{logical indicating if the transposable data were centered.}
}
\author{
Anestis Touloumis
}
\references{
Touloumis, A., Marioni, J.C. and Tavare, S. (2013). Hypothesis Testing for the Covariance Matrix in High-Dimensional Transposable Data with Kronecker Product Dependence Structure. \url{http://arxiv.org/abs/1404.7684}.
}
\examples{
data(VEGFmouse)
## Hypothesis tests for the covariance matrix of the genes (rows).
genestest <- covmat.ts(VEGFmouse,40,"rows",FALSE)
genestest
## Hypothesis tests for the covariance matrix of the tissues (columns).
tissuestest <- covmat.ts(VEGFmouse,40,"columns",FALSE)
tissuestest
}