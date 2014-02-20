\name{trio.power}
\alias{trio.power}
\alias{print.trio.power}
\title{Power and sample size calculation
}
\description{
Computes power for genotypic TDT, allelic TDT or Score test given n trios or required sample size to gain given power.
}
\usage{
trio.power(maf = 0.5, RR = 1.5, alpha = 5*10^(-8), n = NULL, beta = NULL, 
  model = c("additive", "dominant", "recessive"), test = c("gTDT", "Score", "aTDT"))
  
\method{print}{trio.power}(x,digits=4,...)
}
\arguments{
  \item{maf}{
a numeric vector of population frequencies of a mutant allele.
}
  \item{RR}{
a numeric vector of the assumed relative risks for an individual getting a disease with 1 (in case of recessive model 2) mutant alleles compared to the risk of individuals carnying 0 mutant alleles.
}
  \item{alpha}{
a numeric vector of significance levels (Type I Error probability).
}
  \item{n}{
a numeric vector containing number of trios in a study. Must be filled for power calculation. Must not be NULL for sample size calculation.
}
  \item{beta}{
the desired power of the test. Must be filled for power calculation. Must not be NULL for sample size calculation.
}
  \item{model}{
a character containing the genotypic model assumed. Possible values are \code{"additive"}, \code{"dominant"} and \code{"recessive"}. In case of \code{test="aTDT"},
  the standard multiplicative model will be considered. Abbreviations are allowed. Thus, e.g., \code{model = "dom"} will fit a dominant model,
  and \code{model = "r"} a recessive model.
}
  \item{test}{
the chosen test. Must be \code{"aTDT"}, \code{"gTDT"} or \code{"Score"}. Abbreviations are allowed. Thus, e.g., \code{test = "g"} will perform a 
   genotypic TDT, and \code{test = "S"} a Score test.
}
\item{x}{an object of class \code{trio.power}.}
\item{digits}{number of digits that should be printed.}
\item{...}{ignored}
}
\details{
Power and sample size calculation is derived on Knapp (1999). The power or the sample size will be calculated for all combinations of p, RR, alpha, test, model and n or beta.
}
\value{
An object of class trio.power containing the following numeric values or vectors, respectively:
\item{model}{the chosen model}
\item{size}{In case of sample size calculation: calculated sample sizes}
\item{beta}{In case of sample size calculation: desired power}
\item{n}{In case of power calculation: given number of trios}
\item{power}{In case of power calculation: calculated power}
\item{alpha}{Type I error}
\item{test}{the chosen test}
\item{RR}{the relative risks assumed}
\item{p}{the assumed allele frequency}
\item{calc}{the type of calculation}

}
\references{
Knapp, M. (1999). A Note on Power Approximations for the Transmission/Disequilibrium Test. \emph{American Journal of Human Genetics}, 64, 1177-1185.

Neumann, C., Taub, M.A., Younkin, S.G., Beaty, T.H., Ruczinski, I., Schwender, H. (2014). Analytic Power and Sample Size Calculation
for the Genotypic Transmission/Disequilibrium Test in Case-Parent Trio Studies. Submitted.
}

\author{
Christoph Neumann
}

\examples{
# The required samples size to reach of power
# of 0.8 when testing SNPs with minor allele
# frequencies of 0.1 and 0.2 with an additive
# or dominant genotypic TDT and score test
# can be determined by 
trio.power(maf = c(0.1, 0.2), beta = 0.8, model = c("add", "dom"))
}

\keyword{design}
