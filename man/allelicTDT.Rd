\name{allelicTDT}
\alias{allelicTDT}
\alias{print.aTDT}

\title{
Allelic TDT
}
\description{
Performs the allelic Transmission/Disequilibrium Test for each SNP contained in a genotype matrix.
}

\usage{
allelicTDT(mat.snp, size = 50, correct = FALSE)

\method{print}{aTDT}(x, top = 5, digits = 4, ...)
}

\arguments{
  \item{mat.snp}{a numeric matrix in which each column represents a SNP. Each column must be
    a numeric vector of length \eqn{3 * t} representing a SNP genotyped at \eqn{t} trios. Each of the \eqn{t}
    blocks must consist of the genotypes of father, mother, and offspring
    (in this order). The genotypes must be coded by 0, 1, and 2. Missing values are allowed and need to be coded by \code{NA}.
    This matrix might be generated from a ped-file by, e.g., employing \code{\link{ped2geno}}.
}
  \item{size}{the number of SNPs considered simultaneously when computing the parameter estimates.
}
  \item{correct}{should the test statistic be continuity corrected? If \code{FALSE}, \eqn{(b-c)^2 / (b+d)} will be used as test
    statistic, where \eqn{b} and \eqn{c} are the off-diagonal elements of the 2x2-table summarizing the transmitted and
    not transmitted alleles from the heterozygous parents. If \code{TRUE}, \eqn{(|b-c| - 1)^2 / (b+d)} will be used
    as test statistic.}
  \item{x}{an object of class \code{aTDT}, i.e. the output of \code{allelicTDT}.
}
  \item{digits}{number of digits that should be printed.}
  \item{top}{number of interactions that should be printed. If \code{top} is less than or equal to zero, set to \code{NA}, or
    larger than the number of SNPs, then the statistics for all SNPs are printed in the order as they were in the genotype matrix
    used as input into \code{colTDT}. Otherwise, the \code{top} interactions with the smallest p-values are printed.
  }
  \item{...}{ignored.}
}
\value{
An object of class \code{aTDT} containing the following numeric vectors:

\item{stat}{values of the test statistic of the allelic TDT,}
\item{pval}{the corresponding p-values.}
}

\references{
Spielman, R.S., McGinnis, R.E., and Ewens, W.J. (1993). Trsnmmission Test for Linkage Disequilibrium:
The Insulin Gene Region and Insulin-Dependent Diabetes Mellitus (IDDM). \emph{American Journal of Human Genetics},
52, 506-516.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{colTDT}}
}

\examples{
# Load the simulated data for the analysis.
data(trio.data)

# Perform an allelic TDT
a.out <- allelicTDT(mat.test)

# By default, the top 5 SNPs are shown.
# Another number of SNPs, e.g., 10, are displayed by
print(a.out, top=10)

# If the results for all SNPs should be shown (or returned), use
print(a.out, top=0)
}

\keyword{htest}
\keyword{models}
\keyword{array}
