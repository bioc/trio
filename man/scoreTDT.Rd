\name{scoreTDT}
\alias{scoreTDT}
\alias{scoreGxE}
\alias{scoreGxG}
\alias{scoreMaxStat}
\alias{print.scoreTDT}
\alias{print.scoreGxE}
\alias{print.maxScoreTrio}

\title{
Score Tests for SNPs, GxE, and GxG Interactions
}
\description{
Performs score tests for all individual SNPs (\code{scoreTDT}), all interactions of each SNP with an environmental variable (\code{scoreGxE}), 
or all interactions of two SNPs (\code{scoreGxG}) comprised by an input matrix based on the same log-likelihood considered in the corresponding
genotypic TDT, where in \code{scoreGxG} the conditional logistic regression model including only one parameter (for the interaction effect)
is used. 

Additionally, the maximum over the score statistics for testing an additive, dominant, and recessive effect can be
determined using \code{scoreMaxStat}.
}
\usage{
scoreTDT(mat.snp, model = c("additive", "dominant", "recessive"), size = 20)

scoreGxE(mat.snp, env, model = c("additive", "dominant", "recessive"), size = 20, 
   famid = NULL)
   
scoreGxG(mat.snp, model = c("additive", "dominant", "recessive"), genes = NULL, 
   size = 20)
   
scoreMaxStat(mat.snp, size = 20)

\method{print}{scoreTDT}(x, top = 5, digits = 4, ...)

\method{print}{scoreGxE}(x, top = 5, digits = 4, onlyGxE = FALSE, ...)

\method{print}{maxScoreTrio}(x, top = 5, digits = 4, ...)
}

\arguments{
  \item{mat.snp}{a numeric matrix in which each column represents a SNP. Each column must be
    a numeric vector of length \eqn{3 * t} representing a SNP genotyped at \eqn{t} trios. Each of the \eqn{t}
    blocks must consist of the genotypes of father, mother, and offspring
    (in this order). The genotypes must be coded by 0, 1, and 2. Missing values are allowed and need to be coded by \code{NA}.
    This matrix might be generated from a ped-file by, e.g., employing \code{\link{ped2geno}}.
}
  \item{model}{type of model that should be fitted. Abbreviations are allowed. Thus, e.g., \code{model = "dom"} will
    fit a dominant model, and \code{model = "r"} an recessive model. 
}
  \item{size}{the number of models considered simultaneously when computing the parameter estimates.
}
  \item{env}{a vector of length \eqn{t} (see \code{mat.snp}) containing for each offspring the value of a binary environmental
    variable, which must take the values 0 and 1.
}
  \item{famid}{a vector of the same length as \code{env} specifying the family IDs for the corresponding values of the environmental
    variable in \code{env}. Can be used to reorder the vector \code{env} when the order of the trios differs between \code{env} and
    \code{mat.snp}.
}
  \item{genes}{a character vector containing the names of the genes (or LD-blocks or other genetic sets of SNPs) to which the SNPs belong. 
    If specified, only the two-way interactions between SNPs from different genes (or LD-blocks or other genetic sets of SNPs) are tested. 
    If \code{NULL}, all two-way interactions between all possible pairs of SNPs are tested.
}
  \item{x}{an object of class \code{scoreTDT}, \code{scoreGxE}, or \code{maxScoreTrio}, i.e. the output of the function \code{scoreTDT} /
    \code{scoreGxG}, \code{scoreGxE}, or \code{scoreMaxStat}, respectively.
}
  \item{digits}{number of digits that should be printed.}
  \item{top}{number of interactions that should be printed. If the number of interactions is smaller than or equal to 
    \code{top}, then the statistics for all interactions are printed in the order of their computation. Otherwise, 
    the \code{top} interactions with the smallest p-values are printed.
  }
  \item{onlyGxE}{logical indicating whether only the statistics for the parameter of the GxE interaction should be printed.
    If \code{FALSE}, the statistics for both parameters in the model are shown.
}
  \item{...}{ignored.}
}

\value{
For \code{scoreTDT} and \code{scoreGxG}, an object of class \code{scoreTDT} containing numeric vectors
\item{score}{the scores for all SNPs or SNP interactions,}
\item{info}{the denominators of the corresponding score statistics},
\item{stat}{the values of the score statistics for all SNPs or SNP interactions},
\item{pval}{the corresponding p-values computed based on a ChiSquare-distribution with 1 degree of freedom.}
}


\author{Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{\code{\link{colTDT}}, \code{\link{colGxE}}, \code{\link{colTDT2way}}}

\examples{
# Load the simulated data.
data(trio.data)

# A score test can be applied to the SNPs in 
# mat.test by
s.out <- scoreTDT(mat.test)

# By default, an additive mode of inheritance is considered.
# Another mode, e.g., the dominant mode can be considered by
sDom.out <- scoreTDT(mat.test, model = "dominant")

# The test statistic of the MAX score test can be computed by
sMax.out <- scoreMaxStat(mat.test)


# The interaction between a binary environmental factor,
# e.g., the gender, and each SNP in mat.test can be tested
# by setting up the vector containing the value of the
# environmental factor for each trio. If we, e.g., assume
# that the children in the first 50 trios are girls
# and in the remaininge 50 trios boys, then this vector
# can be generated by
sex <- rep(0:1, each = 50)

# and the interaction between sex and each SNP in mat.test
# can be tested with a score test by
sgxe.out <- scoreGxE(mat.test, sex)

# The interactions between all pairs of SNPs in mat.test
# can be tested with a score test by
sgxg.out <- scoreGxG(mat.test)
}


\keyword{htest}
\keyword{array}
