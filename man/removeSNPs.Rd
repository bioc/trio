\name{removeSNPs}
\alias{removeSNPs}
\alias{removeTrios}
\alias{orderSNPs}
\alias{colMAFtrio}
\title{
Remove SNPs or Trios
}
\description{
Functions for removing SNPs with a low minor allele frequency or a high percentage of missing values,
for removing trios in which at least one member shows a high percentage of missing values, for ordering the SNPs
by their position in the genome, and for computing the minor allele frequencies of the SNPs based on only the
genotypes of the parents, where each parent is only used once in this computation, even if this person is part
of more than one of the trios.
}
\usage{
removeSNPs(geno, maf = NA, perc.na = NA)

removeTrios(geno, perc.na = 1)

orderSNPs(geno, map, snp = "SNP", orderBy = c("Chr", "Position"))

colMAFtrio(geno)
}

\arguments{
  \item{geno}{a matrix in genotype format, i.e.\ the output of \code{\link{ped2geno}} or \code{\link{read.pedfile}}
    with \code{p2g} set to \code{TRUE}.
}
  \item{maf}{a numeric value. If specified, i.e.\ not \code{NA}, all SNPs with a minor allele frequency less than
    \code{maf} are removed, where \code{maf} can range from 0 and 0.2. If, e.g., \code{maf = 0}, monomorphic SNPs are removed. 
}
  \item{perc.na}{a numeric value between 0 and 1 specifying a cutoff for the percentage of missing values that a SNP or a subject
    is allowed to have. If more than 100 * \code{perc.na}\% of the genotypes of a SNP or a subject is missing, then this SNP
    or the trio to which this subject belong, respectively, is removed \code{geno}.
}
  \item{map}{a data frame containing the chromosome and the position for all the SNPs in \code{geno}.
}
  \item{snp}{a character string giving the (case-sensitive) name of the column of \code{map} containing the SNP IDs used as
    column names in \code{geno}. 
}
  \item{orderBy}{character string of length 2 specifying the (case-sensitive) names of the columns of \code{map} containing
    the chromosomes and the positions of the SNPs in \code{geno}.
}
}
\value{For \code{removeSNPs}, \code{removeTrios}, and \code{orderSNPs}, a reduced or ordered version of \code{geno}.
For \code{colMAFtrio}, a vector containing the minor allele frequencies of the SNPs in \code{geno}.
}
\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\examples{
# Load the simulated data.
data(trio.data)

# All SNPs with a minor allele frequency smaller than 0.1
# can be removed from mat.test by
mat2 <- removeSNPs(mat.test, maf = 0.1)

# The minor allele frequencies for all SNPs can be 
# determined (based on the genotypes of the parents) by
maf <- colMAFtrio(mat.test)
}

\keyword{manip}

