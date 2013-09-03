\name{getLD}
\alias{getLD}
\alias{getLDlarge}

\title{
 Computation of LD Measures
}
\description{
 While \code{getLD} computes the value of D' and r^2 for each pair of SNPs in a matrix,
 \code{getLDlarge} determines D' and r^2 between each SNP and a user-specified number of SNPs
 closest to the SNP on the corresponding chromosome. Thus, \code{getLDlarge} can be applied
 to much more SNPs than \code{getLD}.
}
\usage{
  getLD(x, which = c("both", "rSquare", "Dprime"), parentsOnly = FALSE,
     iter = 50, snp.in.col = TRUE, asMatrix = FALSE, addVarN = FALSE)
     
  getLDlarge(x, neighbors=25, which=c("both", "rSquare", "Dprime"), 
        parentsOnly=FALSE, iter=50, snp.in.col=TRUE, addVarN=FALSE)     
}

\arguments{
  \item{x}{a numeric matrix consisting of 0, 1, and 2, where it is assumed that
     the values represent the numbers of minor alleles that the SNPs show. Missing values
     are allowed. By default, each column represents a SNP and each row a subject.
     This can be changed by setting \code{snp.in.col = FALSE}. It is assumed that the SNPs are ordered
     by their position on the considered chromosome.
}
  \item{neighbors}{positive integer specifying the number of neighbors of a SNP (in both directions) on a chromosome
    for which D' or r^2 should be computed. Thus, for each SNP (except for the SNPs in the first and last \code{neighbors}
    columns of \code{x}), 2 * \code{neighbors} r^2 or D' values are computed.
}
  \item{which}{which LD measures should be computed? Either \code{"rSquare"}, or \code{"Dprime"},
    or the values of \code{"both"} measures are computed. The latter is the default.
}
  \item{parentsOnly}{logical indicating whether only the genotypes of the parents, i.e.\ rows 1, 2,
    4, 5, ... of \code{x}, should be used in the computation of the LD measures when \code{x}
    is in genotype format and contains case-parent trio data (see \code{\link{ped2geno}} and \code{\link{read.pedfile}}).
    If \code{FALSE} (default), all rows are used in the determination of the pairwise LD measure.
} 
  \item{iter}{integer specifying how many iterations are used in the procedure of Hill (1974)
    which is used to estimate D.
}
  \item{snp.in.col}{logical indicating whether each column of \code{x} represents a SNP (and each row a
    subject). If \code{FALSE}, each row represents a SNP (and each column a subject).
}
  \item{asMatrix}{logical indicating whether the LD values are returned as a \eqn{m} x \eqn{m} matrix,
    where \eqn{m} is the number of SNPs. If \code{FALSE}, the LD values are returned as a vector
    of length \eqn{m * (m - 1) / 2}.
}
  \item{addVarN}{logical indicating whether for each pair of SNPs the number of non-missing values
    and the variance estimates of D' proposed by Zabaleta et al. (1997) should be added to the output. The
    variance estimates are required for the identification of LD-blocks with \code{\link{findLDblocks}}.
}
}

\value{
An object of class \code{getLD} or \code{getLDlarge} consisting (depending of the specification of \code{which}) the
D' (\code{Dprime}) or r^2 (\code{rSquare}) values for each SNP pair, and (depending of the specification
of \code{addVarN}) the variance estimates for D' (\code{varDprime}) and the numbers of non-missing values
(\code{n}). Furthermore, the names of the SNPs (\code{rn}) will be added (in \code{getLD}, if \code{asMatrix = FALSE}).

}
\references{
  Hill, W.O. (1974). Estimation of Linkage Disequilibrium in Randomly Mating Populations. 
    \emph{Heredity}, 33, 229-239.
  
  Zapata, C., Alvarez, G., and Carollo, C. (1997). Approximate Variance of the Standardized Measure of Gametic
     Disequilibrium D'. \emph{American Journal of Human Genetics}, 61, 771-774.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{plot.getLD}}, \code{\link{findLDblocks}}
}

\examples{
# Load the simulated data.
data(trio.data)

# The values of Dprime and Rsquare for each pair of SNPs 
# in LDdata can be computed by
ld.out <- getLD(LDdata)

# By default, the LD measures are returned as a vector.
# If they should be returned as a matrix, then use
ld.out2 <- getLD(LDdata, asMatrix = TRUE)
}

\keyword{array}
\keyword{univar}
