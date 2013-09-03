\name{colTDTsam}
\alias{colTDTsam}
\alias{colTDTebam}
\alias{gtdt.stat}
\alias{gtdt.ebam}

\title{
SAM and EBAM for Trio Data
}
\description{
Performs a Significance Analysis of Microarrays (SAM; Tusher et al., 2001) or an Empirical Bayes Analysis of Microarrays (EBAM;
Efron et al., 2001), respectively, based on the genotypic transmission/disequilibrium test statistic.
}
\usage{
colTDTsam(mat.snp, model = c("additive", "dominant", "recessive", "max"), 
   approx = NULL, B = 1000, size = 10, chunk = 100, rand = NA)
   
colTDTebam(mat.snp, model = c("additive", "dominant", "recessive", "max"), 
   approx = NULL, B = 1000, size = 10, chunk = 100, 
   n.interval = NULL, df.ratio = 3, df.dens = 3, knots.mode = TRUE, 
   type.nclass = c("wand", "FD", "scott"), fast = FALSE, rand = NA)
}

\arguments{
  \item{mat.snp}{a matrix in genotype format, i.e. a numeric matrix in which each column is
    a vector of length \eqn{3 * t} representing a SNP genotyped at \eqn{t} trios. Each of the \eqn{t}
    blocks of rows in \code{mat.snp} must consist of the genotypes of father, mother, and offspring
    (in this order), where the genotypes must be coded by 0, 1, and 2. Missing values are allowed and need to be coded by \code{NA}.
    This matrix might be generated from a data frame in ped format by, e.g., employing \code{\link{ped2geno}}.
}
  \item{model}{type of genetic mode of inheritance that should be considered. Either \code{"additive"} (default), \code{"dominant"},
    \code{"recessive"}, or \code{"max"}. If \code{model = "max"}, the maximum over the gTDT statistics for testing an additive, dominant,
    and recessive model is used as gTDT statistic. Abbreviations are allowed. Thus, e.g., \code{model = "dom"} will
    fit a dominant model, and \code{model = "r"} an recessive model. 
}
  \item{approx}{logical specifying whether the null distribution should be approximated by a \eqn{\chi^2}{ChiSquare}-distribution
     with one degree of frredom. If \code{approx = FALSE}, the null distribution is estimated based on a permutation method.
     If not specified, i.e. \code{NULL}, \code{approx} is set to \code{TRUE}, when an additive, dominant,
     or recessive mode of inheritance is considered, and \code{approx = FALSE}, when \code{model = "max"}.
     If \code{model = "max"}, it is not allowed to set \code{approx = TRUE}.
}
  \item{B}{number of permutations used in the estimation of the null distribution, and thus, the computation of the null statistics.
     Ignored if \code{approx = TRUE}.
}
  \item{size}{number of SNPs considered simultaneously when computing the gTDT statistics.
}
  \item{chunk}{number of permutations considered simultaneously in the permutation procedure.
}
  \item{n.interval}{the number of intervals used in the logistic regression with
     repeated observations for estimating the ratio of the null density to the density of the observed
     gTDT values in an EBAM analysis (if \code{approx = FALSE}), or in the Poisson regression used to estimate
     the density of the observed gTDT values (if \code{approx = TRUE}). For details, see Efron et al., 2001,
     or Schwender and Ickstadt, 2008, respectively.
     If \code{NULL}, \code{n.interval} is determined by the maximum of 139 (see Efron et al., 2001) and the 
     number of intervals estimated by the method specified by \code{type.nclass}.
}
  \item{df.ratio}{integer specifying the degrees of freedom of the natural cubic
     spline used in the logistic regression with repeated observations for estimating the ratio of the null
     density to the density of the observed gTDT values in an EBAM analysis. Only used when \code{approx} is
     set to \code{FALSE}. 
}
  \item{df.dens}{integer specifying the degrees of freedom of the natural cubic
     spline used in the Poisson regression to estimate the density of the observed gTDT values in an EBAM analysis.
     Only used when \code{approx} is set to \code{TRUE}.
}
  \item{knots.mode}{logical specifying whether the \code{df.dens} - 1 knots of the natural cubic spline are centered around the
     mode and not the median of the density when fitting the Poisson regression model to estimate
     the density of the observed gTDT values in an EBAM analysis. Only used when \code{approx} is set to \code{TRUE}.
     For details on this density estimation, see \code{\link[siggenes]{denspr}}. 
}
  \item{type.nclass}{character string specifying the procedure used to estimate the
     number of intervals of the histogram used in the logistic regression with repeated observations or the Poisson regression,
     respectively (see \code{n.interval}).  Can be either \code{"wand"} (default), \code{"FD"}, or \code{"scott"}. 
     Ignored if \code{n.interval} is specified. For details, see \code{\link[siggenes]{denspr}}.
}
  \item{fast}{logical specifying whether a crude estimate for the number of permuted test scores larger than the respective
     observed gTDT value should be used. If \code{FALSE}, the exact number of permuted test scores larger than the respective
     observed gTDT value is computed. 
}
  \item{rand}{numeric value. If specified, i.e. not \code{NA}, the random number generator
     will be set into a reproducible state.
}
}

\value{
The output of \code{colTDTsam} or \code{colTDTebam} is an object of class \code{SAM} or \code{EBAM}, respectively. All the
features implemented in the \code{R} package \code{siggenes} for an SAM or EBAM analysis, respectively, can therefore be
used in the SAM or EBAM analysis of case-parent trio data implemented in \code{colTDTsam} or \code{colTDTebam}, respectively.
For details, see \code{\link[siggenes]{sam}} or \code{\link[siggenes]{ebam}}, respectively.
}

\references{
Efron, B., Tibshirani, R., Storey, J.D., and Tusher, V. (2001). 
   Empirical Bayes Analysis of a Microarray Experiment, \emph{Journal of the American Statistical Association}, 
   96, 1151-1160.

Schwender, H. and Ickstadt, K. (2008). Empirical Bayes Analysis of Single Nucleotide Polymorphisms. 
   \emph{BMC Bioinformatics}, 9, 144.

Schwender, H., Taub, M.A., Beaty, T.H., Marazita, M.L., and Ruczinski, I. (2011). 
Rapid Testing of SNPs and Gene-Environment Interactions in Case-Parent Trio Data Based on 
Exact Analytic Parameter Estimation. \emph{Biometrics}, 68, 766-773.

Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance Analysis of Microarrays
   Applied to the Ionizing Radiation Response. \emph{Proceedings of the National Academy of Science of the 
   United States of America}, 98, 5116-5121.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{colTDT}}, \code{\link{colTDTmaxStat}}, \code{\link[siggenes]{sam}}, \code{\link[siggenes]{ebam}},
\code{\link[siggenes]{SAM-class}}, \code{\link[siggenes]{EBAM-class}}
}

\examples{
# Load the simulated data.
data(trio.data)

# Perform a Significance Analysis of Microarrays (SAM).
sam.out <- colTDTsam(mat.test)

# By default an additive mode of inheritance is considered.
# If another mode, e.g., the dominant mode, should be 
# considered, then this can be done by
samDom.out <- colTDTsam(mat.test, model="dominant")

# Analogously, an Empirical Bayes Analysis of Microarrays based
# on the genotypic TDT can be performed by
ebam.out <- colTDTebam(mat.test)
}

\keyword{htest}
\keyword{array}
\keyword{models}
