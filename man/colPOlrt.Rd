\name{colPOlrt}
\alias{colPOlrt}
\alias{colTAT}
\alias{print.polrt}
\alias{print.tat}
\title{
Parent-of-Origin Tests
}
\description{
Computes the test statistics and the corresponding p-values either for the Parent-of-Origin Likelihood Ratio Test proposed by Weinberg (1999)
or the Transmission Asymmetry Test proposed by Weinberg et al. (1998).
}

\usage{
colPOlrt(mat.snp, size = 20)
   
colTAT(mat.snp, stratified = FALSE, size = 50, bothHet = 0) 

\method{print}{polrt}(x, top = 5, digits = 4, ...)

\method{print}{tat}(x, top = 5, digits = 4, ...)
}

\arguments{
  \item{mat.snp}{a numeric matrix in which each column represents a SNP. Each column must be
    a numeric vector of length \eqn{3 * t} representing a SNP genotyped at \eqn{t} trios. Each of the \eqn{t}
    blocks must consist of the genotypes of father, mother, and offspring
    (in this order). The genotypes must be coded by 0, 1, and 2. Missing values are allowed and need to be coded by \code{NA}.
    This matrix might be generated from a ped-file by, e.g., employing \code{\link{ped2geno}}.
}
  \item{size}{the number of SNPs considered simultaneously when computing the test statistics.
}
  \item{stratified}{should also test statistics and p-values stratified by paternal and maternal transmission be computed?
}
  \item{bothHet}{a numeric value between 0 and 1 specifying how trios in which both parents are heterozygous are weighted in
     determination of the TAT statistic. By default, such trios are ignored (as proposed by Weinberg, 1999). If \code{bothHet = 1},
     such trios are treated in the same way as trios with one heterozygous parent. Other values (e.g., \code{bothHet = 0.5})
     are also sometimes used for \code{bothHet}.
}
  \item{x}{an object of class \code{polrt} or \code{tat}, i.e. the output of \code{colPOlrt} or \code{colTAT}, respectively.
}
  \item{digits}{number of digits that should be printed.}
  \item{top}{number of interactions that should be printed. If \code{top} is less than or equal to zero, set to \code{NA}, or
    larger than the number of SNPs, then the statistics for all SNPs are printed in the order as they were in the genotype matrix
    used as input into \code{colTDT}. Otherwise, the \code{top} interactions with the smallest p-values are printed.
  }
  \item{...}{ignored.}}

\value{
For \code{colPOlrt}, an object of class \code{polrt} consisting of the following numeric vectors:

\item{stat}{the values of the test statistic of the likelihood ratio test for all SNPs in \code{mat.snp},}
\item{pval}{the corresponding p-values,}
\item{full}{the values of the maximized likelihoods of the full models containing also a parameter for the parent-of-origin effect,}
\item{red}{the values of the maximied likelihoods of the reduced models not containing this parameter.}

For \code{colTAT}, an object of class \code{tat} consisting of the following numeric vectors:

\item{stat}{the values of the test statistic of transmission asymmetry test for all SNPs in \code{mat.snp},}
\item{pval}{the corresponding p-values,}
\item{usedTrios}{the number of trios affecting the determination of the TAT statistic,}

and if \code{stratified = TRUE}

\item{matStrat}{a matrix with four columns containing the number of minor alleles transmitted and not-transmitted by heterozygous 
fathers and mothers,}
\item{statPaternal}{a numeric vector containing for each SNP the value of the test statistic for testing whether 
the numbers of paternal transmissions and non-transmissions differ,}
\item{pvalPaternal}{the p-values corresponding to \code{statPaternal},} 
\item{statMaternal}{a numeric vector containing for each SNP the value of the test statistic for testing whether 
the numbers of maternal transmissions and non-transmissions differ,}
\item{pvalMaternal}{the p-values corresponding to \code{statMaternal}.} 
}

\references{
Weinberg, C.R.,Wilcox, A.J., and Lie, R.T. (1998). A Log-Linear Approach to Case-Parent-??Triad Data: 
Assessing Effects of Disease Genes that act Either Directly or Through Maternal Effects and that
may be Subject to Parental Imprinting. \emph{American Journal of Human Genetics}, 62, 969-978.

Weinberg, C.R. (1999). Methods for Detection of Parent-of-Origin Effects in Genetic Studies of
Case-Parents Triads. \emph{American Journal of Human Genetics}, 65, 229-235.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{colTDT}}, \code{\link{ped2geno}}
}

\examples{
# Load the simulated data.
data(trio.data)

# The Parent-of-Origin Likelihood Ratio Test can be applied
# to the SNPs in mat.test by
po.out <- colPOlrt(mat.test)

# The Transmission Asymmetry Test can be applid to the SNPs
# in mat.test by 
tat.out <- colTAT(mat.test)

# By default, statistics for the top 5 SNPs are displayed.
# If another number of SNPs, say 10, should be displayed,
# then this can be done by
print(po.out, top = 10)

# The statistics for all SNPs (not ordered by their 
# significance) can be obtained by
print(po.out, top = 0)
}


\keyword{htest}
\keyword{models}
\keyword{array}
