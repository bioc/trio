\name{ntrios2Dom}
\alias{ntrios2Dom}
\alias{ntrios2Rec}

\title{
Frequency of Dominant or Recessive Genotype Combinations
}
\description{
Transform the numbers of trios showing the different combinations of genotypes in trios coded by 0, 1, and 2 (i.e. considering an additive mode of inheritance) to the numbers of trios showing the different genotype combinations when the three genotypes of a SNP are coded for a dominant or a recessive mode of inheritance.
}

\usage{
ntrios2Dom(matNumber, check = TRUE, quiet = FALSE)

ntrios2Rec(matNumber, check = TRUE, quiet = FALSE)
}


\arguments{
  \item{matNumber}{a numeric matrix or list of two matrices containing the numbers of trios showing the different combinations of genotypes of the offsprings and their parents, i.e. the output of \code{\link{colNtrios}}. Each column of each matrix corresponds to a genotype combination, and each row to a SNP.
}
  \item{check}{logical specifying whether it should be checked if it is plausible that \code{matNumber} is the output of \code{\link{colNtrios}}. Meant for internal purposes only. So the default should actually not be changed.
}
  \item{quiet}{logical indicating whether a message should be surpressed when only the seven numbers of trios contributing to the likelihood estimation considered in the genotypic TDT have been computed.
}
}




\value{
If in \code{\link{colNtrios}} the argument \code{env} has not been specified, a matrix with the numbers of trios, where each column of this matrix corresponds to a genotype combination and each row to a SNP. Depending on whether \code{ntrios2Dom} or \code{ntrios2Rec} has been applied, the column names of this matrix start with a D or R (for dominant or recessive mode of inheritance), respectively, followed by the coding for the offspring and then the three pseudo controls (see Table 2 or 3, respectively, in Schwender et al., 2011).

If \code{env} has been specified in \code{\link{colNtrios}}, a list consisting of two matrices as the one described above is returend, one matrix for \code{env = 0}, and the other for \code{env = 1}.
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
\code{\link{colNtrios}}, \code{\link{colTDT}}
}

\examples{
# Load the simulated data.
data(trio.data)

# Compute the numbers of trios.
mat.ntrios <- colNtrios(mat.test) 

# Transform these numbers to the numbers
# when considering genotypes coding for a
# dominant mode of inheritance.
ntrios2Dom(mat.ntrios)
}

\keyword{array}
