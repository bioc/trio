\name{colGxGPerms}
\alias{colGxGPerms}
\alias{compPermTDT2way}
\title{
Permutation-Based gTDT for Two-Way Interactions
}
\description{
Computes the original and permuted values of the test statistic of the gTDT test as proposed by Cordell (2002) for each
interaction between the pairs of SNPs in mat.snp.
}
\usage{
colGxGPerms(mat.snp, n.perm = 1000, genes = NULL, col.out = NULL, 
   warnError = TRUE, verbose = TRUE, rand = NA)
}
\arguments{
  \item{mat.snp}{a numeric matrix in which each column represents a SNP. Each column must be
    a numeric vector of length \eqn{3 * t} representing a SNP genotyped at \eqn{t} trios. Each of the \eqn{t}
    blocks must consist of the genotypes of father, mother, and offspring
    (in this order). The genotypes must be coded by 0, 1, and 2. Missing values are allowed and need to be coded by \code{NA}.
    This matrix might be generated from a ped-file by, e.g., employing \code{\link{ped2geno}}.
}
  \item{n.perm}{number of permutations of the response for which the permuted values of the test statistic should be computed.
}
  \item{genes}{a character vector containing the names of the genes to which the SNPs belong. If specified, only the two-way interactions
    between SNPs from different genes are tested. If \code{NULL}, all two-way interactions between all possible pairs of SNPs are tested.
}
  \item{col.out}{the output of \code{\link{colGxG}} with \code{epistatic = TRUE} (which is the default in \code{\link{colGxG}}).
    If \code{NULL}, \code{compPermTDT2way} computes the values of the test statistic for the original permutation of the response.
}
  \item{warnError}{logical indicating whether the statistics for the gTDT should be returned as \code{NA} if the fitting of the
    conditional logistic regression model fails. This might in particular happen when the two considered SNPs are in (strong) LD.
}
  \item{verbose}{logical indicating whether some information on what is currently computed should be printed.
}
  \item{rand}{numeric value. If specified, the random number generator is set into a reproducible state.
}
}
\value{
A list consisting of
\item{stat}{a numeric vector containing the original values of the test statistic,}
\item{permStat}{a numeric matrix containing the permuted values of the test statistic,}
\item{y.perm}{a matrix containing the permutations of the response.}
}

\references{
Cordell, H. J. (2002). Epistasis: What it Means, what it Doesn't mean, and Statistical Methods to Detect it in Humans.
Human Molecular Genetics, 11, 2463-2468.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{colGxG}}
}

\examples{
# Load the simulated data.
data(trio.data)

# Cordell's LRT for all pairs of SNPs in mat.test can be performed
# and the values of the LRT statistic for 10 permutations of the
# case-pseudo-controls status can be computed by
gxg <- colGxGPerms(mat.test, n.perm = 10)

# where we here consider only 10 permutations to keep the computing
# time of this example small. Usually, at least a few thousand 
# permutations should be considered.
}


\keyword{htest}
\keyword{models}
\keyword{array}
