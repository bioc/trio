\name{colEMlrt}
\alias{colEMlrt}
\alias{print.colEMlrt}

\title{
EM Likelihood Ratio Test
}
\description{
Performs the Expectation Maximimization Likelihood Ratio Test (EMLRT) proposed by Weinberg (1999) for each SNP in a 
matrix in genotype format. 
}
\usage{
colEMlrt(mat.snp, model = c("general", "dominant", "recessive"), maternal = FALSE, 
     parentMissing = c("father", "mother", "either"), iter = 40, eps = 10^-16)

\method{print}{colEMlrt}(x, top = 5, digits = 4, ...)
}

\arguments{
  \item{mat.snp}{a numeric matrix in which each column represents a SNP. Each column must be
    a numeric vector of length \eqn{3 * t} representing a SNP genotyped at \eqn{t} trios. Each of the \eqn{t}
    blocks must consist of the genotypes of father, mother, and offspring
    (in this order). The genotypes must be coded by 0, 1, and 2. Missing values are allowed and need to be coded by \code{NA}.
    This matrix in genotype format might be generated from a ped-file by, e.g., employing \code{\link{ped2geno}}.
}
  \item{model}{a character string specifying the genetic effect that should be considered in the Poisson regression. 
	   By default, the general model proposed by Weinberg (1999) is fitted. 
     Alternatively, a dominant or recessive mode of inheritance might be considered by setting \code{model = "dominant"}
		 or \code{model = "recessive"}, respectively. Abbreviations such as \code{"dom"} or \code{"rec"}
		are also accepted.
}
  \item{maternal}{logical specifying whether parameters for a maternal effects should be added to the Poisson regression model. 
	   If \code{TRUE} and \code{model = "general"}, the most general model described by Weinberg (1999)
		 is used.
}
\item{parentMissing}{a character string specifying whether the genotype of the \code{father}, the \code{mother}, or \code{either} 
   the mother or the father is allowed to be missing. By default, only the genotype of the father is allowed to be missing in a trio.
	 Thus, in this case, all complete trios and all trios in which data are available for the father and the offspring are used in the
	 testing of the considered SNP. If \code{parentMissing = "either"}, all trios in which the genotype of the offspring or the genotypes
	 of both parents are missing for a particular SNP are ignored in the analysis of this SNP.
}
  \item{iter}{a non-negative numeric value specifying the maximum number of iterations that should be used in the application of the expectation maximization algorithm.
}
  \item{eps}{a non-negative small numeric value specifying the accuracy for the stopping criterion of the EM algorithm.
	   If the sum of the squared differences of the estimated expected numbers of trios over the 15 possible genotype combinations
		in trios between two consecutive iterations of the EM algorithm is smaller than \code{eps}, the EM algorithm stops. Often, less than ten
		 genotype combinations are required.
}
  \item{x}{an object of class \code{colEMlrt}, i.e. the output of \code{colEMlrt}.
}
  \item{digits}{number of digits that should be printed.}
  \item{top}{number of SNPs that should be printed. If \code{top} is less than or equal to zero or
    larger than the total number of SNPs in \code{mat.snp}, then the statistics for all SNPs are printed in the order as they were in the genotype matrix
    used as input into \code{colEMlrt}. Otherwise, the \code{top} interactions with the smallest p-values are printed.
  }
  \item{...}{ignored.}
}
\details{
While in functions such as \code{\link{colTDT}} all trios in which the genotype of one or more of the members of this trio 
is missing for a particular SNP are removed from the analysis of this SNP, the procedure proposed by Weinberg (1999) can handle missing 
genotypes by employing an expectation maximization (EM) algorithm to estimate the expected numbers of
trios for the 15 different genotype combinations
possible in a trio (when considering the genotypes of mothers and fathers individually) and a likelihood ratio test based
on two nested Poisson regression models using the estimated expected numbers of trios as outcome.

}

\value{
An object of class \code{colEMlrt} consisting of the following numeric vectors:

\item{stat}{the values of the test statistic of the likelihood ratio test for all SNPs in \code{mat.snp},}
\item{pval}{the corresponding p-values,}
\item{ll.full}{the values of the maximized likelihoods of the full models,}
\item{ll,red}{the values of the maximied likelihoods of the reduced models not containing the
   parameter(s) corresponding to the genetic \code{model} of interest.}
}
\references{
Weinberg, C.R. (1999). Allowing for Missing Parents in Genetic Studies of Case-Parent Triads. 
\emph{American Journal of Human Genetics}, 64, 1186-1193.
}
\author{
Philipp Berger, \email{philipp.berger@hhu.de}
}

\seealso{
\code{\link{colTDT}}, \code{\link{ped2geno}}
}

\examples{
# Load the simulated data.
data(trio.data)

# The EM Likelihood Ratio Test can be applied
# to the SNPs in mat.test by
em.out <- colEMlrt(mat.test)

# If a dominant mode of inheritance is of interest,
# the corresponding EM Likelihood Ratio Test can be
# performed by 
emd.out <- colEMlrt(mat.test, model="dominant")

# By default, statistics for the top 5 SNPs are displayed.
# If another number of SNPs, say 10, should be displayed,
# then this can be done by
print(em.out, top = 10)

# The statistics for all SNPs (not ordered by their 
# significance) can be obtained by
print(em.out, top = 0)
}

\keyword{htest}
\keyword{models}
\keyword{array}
