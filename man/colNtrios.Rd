\name{colNtrios}
\alias{colNtrios}

\title{
Frequency of Genotype Combinations
}
\description{
Computes the numbers of trios showing the ten different genotype combinations for each SNP in a matrix in genotype format.
}

\usage{
colNtrios(mat.snp, env = NULL, onlyContributing = FALSE, famid = NULL,
   size = 50)
}


\arguments{
  \item{mat.snp}{a numeric matrix in which each column represents a SNP. Each column must be a numeric vector of length \eqn{3 * t} representing a SNP genotyped at \eqn{t} trios. Each of the \eqn{t}
    blocks must consist of the genotypes of father, mother, and offspring (in this order). The genotypes must be coded by 0, 1, and 2. Missing values are allowed and need to be coded by \code{NA}.
    This matrix in genotype format might be generated from a ped file by, e.g., employing \code{\link{ped2geno}}.
}
  \item{env}{a vector of length \eqn{t} (see \code{mat.snp}) containing for each offspring the value of a binary environmental variable, which must take the values 0 and 1. If specified the numbers of trios showing the ten different genotype combinations are computed separately for the trios for which \code{env = 0} and for the trios for which \code{env = 1}.
}
  \item{onlyContributing}{if \code{TRUE}, only the seven numbers of trios that contribute to the likelihood estimation in the performance of the genotypic TDT, e.g., with \code{\link{colTDT}}, are determined. If \code{FALSE} (which is the default), all ten numbers of trios for the ten different genotype combinations are computed.}
	  \item{famid}{a vector of the same length as \code{env} specifying the family IDs for the corresponding values of the environmental variable in \code{env}. Can be used to reorder the vector \code{env} when the order of the trios differs between \code{env} and \code{mat.snp}. Ignored if \code{env} is not specified.
}
  \item{size}{the number of SNPs considered simultaneously when computing the parameter estimates. Ignored if \code{fast = FALSE}.
}
}


\details{
When considering SNPs coded by the numbers of minor alleles so that the genotypes of the SNPs are coded by 0, 1, and 2, there exist ten possible combinations of the genotypes of an offspring and its parents that this trio can show (see Table 1 in Schwender et al., 2011). Seven of these genotype combinations contribute to the likelihood estimation considered when performing a genotypic TDT.

Depending on the specification of the argument \code{onlyContributing}, \code{colNtrios} computes for each SNP in \code{mat.snp} the numbers of trios showing the ten or seven different genotype combinations, respectively.

In \code{colNtrios}, an additive mode of inheritance is considered. If the numbers of trios showing the different genotype combinations that exist when considering a dominant or recessive mode of inheritance are of interest (see, e.g., Table 2 or 3, respectively, in Schwender et al., 2011), then \code{\link{ntrios2Dom}} or \code{\link{ntrios2Rec}}, respectively, can be applied to the output of \code{colNtrios}.
} 

\value{
If \code{env} is not specified, a matrix with the numbers of trios for the different genotype combinations, where the names of the columns corresponding to the different genotype combinations are given	by an starting "G" for genotype followed by the genotypes of the parents and the offspring (as in Table 1 of Schwender et al., 2011).

If \code{env} is specified, a list consisting of two matrices as the one described above, one matrix for \code{env = 0}, and the other for \code{env = 1}.
}

\references{
Schwender, H., Taub, M.A., Beaty, T.H., Marazita, M.L., and Ruczinski, I. (2011). 
Rapid Testing of SNPs and Gene-Environment Interactions in Case-Parent Trio Data Based on 
Exact Analytic Parameter Estimation. \emph{Biometrics}, 68, 766-773. 
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{colTDT}}, \code{\link{ntrios2Dom}}
}

\examples{
# Load the simulated data.
data(trio.data)

# Compute the numbers of trios.
mat.ntrios <- colNtrios(mat.test) 

# A genotypic TDT can then be applied for each SNP 
# in mat.test by
colTDT(matNumber = mat.ntrios)

# This leads to the same results as
colTDT(mat.test)

# If env is also specified, e.g., by
facEnv <- rep(0:1, each = 50)
list.ntriosEnv <- colNtrios(mat.test, env = facEnv)

# then a genotypic TDT for gene-environment interactions
# can be performed by
gxe.Nout <- colGxE(listNumber = list.ntriosEnv, 
    env = facEnv)
}

\keyword{models}
\keyword{array}
