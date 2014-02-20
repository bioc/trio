\name{print.colGxE}
\alias{print.colGxE}
\alias{getGxEstats}
\title{
Printing and Storing of colGxE objects
}
\description{
\code{Print}s the statistics computed with \code{\link{colGxE}}. \code{getGxEstats} generates a data frame
containing these statistics.
}
\usage{
\method{print}{colGxE}(x, top = 5, digits = 4, onlyGxE = FALSE, ...)

getGxEstats(x, top = NA, sortBy = c("none", "gxe", "lrt2df", "wald2df", "lrt1df", "g"))
}

\arguments{
  \item{x}{an object of class \code{colGxE}, i.e. the output of the function \code{\link{colGxE}}.
}
  \item{top}{number of top interactions that should be printed or stored in a data frame. If \code{top} is set to \code{NA}, \code{0}, 
    or to a value that is negative of larger than the number of interactions, then the statistics for all interactions are printed or stored
    in the same order as they were in the genotype matrix \code{mat.snp} used in \code{\link{colGxE}}. Otherwise, 
    the \code{top} interactions with the smallest p-values are printed or stored, where \code{print} uses the p-values of the GxE effect
    to order the interactions, while in \code{generateGxEstats} the p-values of test specified by \code{sortBy} are employed. Ignored
    if \code{sortBy = "none"}.
}
  \item{onlyGxE}{logical indicating whether only the statistics for the parameter of the GxE interaction should be printed.
    If \code{FALSE}, the statistics for both parameters in the model as well as the relative risks for the exposed trios and statistics for
    the 2 df likelihood ratio test and the 2 df Wald test (if these relative risks and statistics were computed by \code{\link{colGxE}})
    are shown.
}
  \item{digits}{number of digits that should be printed.}
  \item{...}{ignored.}
  \item{sortBy}{character string specifying by the p-value of which test the SNPs should be sorted. If \code{"none"} (default), the
    SNPs are not sorted and the SNPs are in the same order as in the genotype matrix used to specify \code{mat.snp} in \code{\link{colGxE}}.}
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
\code{\link{colGxE}}
}

\examples{
# Load the simulated data for the analysis.
data(trio.data)

# Set up a vector with the binary environmental variable.
# Here, we consider the gene-gender interactions and
# assume that the children in the first 50 trios are
# girls, and the remaining 50 are boys.
sex <- rep(0:1, each = 50)

# Test the interaction of sex with each of the SNPs in mat.test
gxe.out <- colGxE(mat.test, sex)

# By default, the statistics are shown for the parameters of
# the top 5 GxE interactions and the parameters of the
# corresponding SNPs. 
gxe.out

# If the top 10 GxE interactions should be displayed, then this
# can be done by
print(gxe.out, top = 10)

# The statististics for all GxE interactions (and SNPs) are
# shown, when calling
print(gxe.out, top = 0)

# If only the statistics for the GxE parameters, but not for
# the SNPs should be displayed, then use
print(gxe.out, onlyGxE = TRUE)

# A convenient way to generate a data frame with all the statistics
# computed by colGxE either for the top SNPs or for all SNPs (here,
# the top 10 SNPs) ordered by the p-values of one of the considered 
# tests, e.g., the 2 df likelihood ratio test, is
dat.top3 <- getGxEstats(gxe.out, top = 10, sortBy = "lrt2df")
}


\keyword{print}
\keyword{utilities}
