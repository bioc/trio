\name{trioFS}
\alias{trioFS}
\alias{trioFS.default}
\alias{trioFS.trioPrepare}
\alias{trioFS.formula}

\title{Trio Feature Selection
}
\description{Performs a trioFS (trio Feature Selection) analysis as proposed by Schwender et al. (2011) based on bagging/subsampling
with base learner trio logic regression (Li et al., 2011).
}
\usage{
\method{trioFS}{default}(x, y, B = 20, nleaves = 5, replace = TRUE, sub.frac = 0.632, 
    control = lrControl(), fast = FALSE, addMatImp = TRUE, addModels = TRUE, 
    verbose = FALSE, rand = NA, ...)

\method{trioFS}{trioPrepare}(x, ...)

\method{trioFS}{formula}(formula, data, recdom = TRUE, ...)
}

\arguments{
  \item{x}{either an object of class \code{trioPrepare}, i.e. the output of \code{\link{trio.prepare}}, or 
      a binary matrix consisting of zeros and ones. If the latter, then each column of \code{x} must correspond to a binary variable 
     (e.g., codng for a dominant or a recessive effect of a SNP), and each row to a case or a pseudo-control,
     where each trio is represented by a block of four consecutive rows of \code{x} containing the data for the case
     and the three matched pseudo-controls (in this order) so that the first four rows of \code{x} comprise
     the data for the first trio, rows 5-8 the data for the second trio, and so on. Missing values are not
     allowed. A convenient way to generate this matrix is to use the function \code{\link{trio.prepare}}. Afterwards, \code{trioLR} 
     can be directly applied to the output of \code{\link{trio.prepare}}.
}
  \item{y}{a numeric vector specifying the case-pseudo-control status for the observations in \code{x} (if \code{x} is a binary matrix).
     Since in trio logic regression, cases are coded by a \code{3} and pseudo-controls by a \code{0},
     \code{y} is given by \code{rep(c(3, 0, 0, 0), n.trios)}, where \code{n.trios} is
     the number of trios for which genotype data is stored in \code{x}. Thus, the length of \code{y}
     must be equal to the number of rows in \code{x}. No missing values are allowed in \code{y}.
     If not specified, \code{y} will be automatically generated. 
}
  \item{B}{number of bootstrap samples or subsamples used in \code{trioFS}
}
  \item{nleaves}{maximum number of leaves, i.e.\ variables, in the logic tree considered in each of the \code{B}
     trio logic regression models (please note in trio logic regression the model consists only of one logic tree).
}
  \item{replace}{should sampling of the trios be done with replacement? If 
     \code{TRUE}, a Bootstrap sample of size \code{n.trios} is drawn
     from the \code{n.trios} trios in each of the \code{B} iterations. If
     \code{FALSE}, \code{ceiling(sub.frac * n.trios)} of the trios
     are drawn without replacement in each iteration.}
  \item{sub.frac}{a proportion specifying the fraction of trios that
     are used in each iteration to fit a trio logic regression model if \code{replace = FALSE}.
     Ignored if \code{replace = TRUE}.}
  \item{control}{a list of control parameters for the search algorithms and the logic trees considered when fitting the
     trio logic regression model, where the parameters for an MC logic regression are ignored. For details and the parameters,
     see \code{\link{lrControl}}, which is the function that should be used to specify \code{control}.}
  \item{fast}{should a greedy search be used instead of simulated annealing, i.e. the standard 
     search algorithm in (trio) logic regression?
}
  \item{addMatImp}{should the matrix containing the improvements due to the interactions
     in each of the iterations be added to the output, where the importance of each interaction
     is computed by the average over the \code{B} improvements due to this interaction?
}
  \item{addModels}{should the \code{B} trio logic regression models be added to the output
}
  \item{verbose}{should some comments on the progress the \code{trioFS} analysis be printed?
}
  \item{rand}{positive integer. If specified, the random number generator is set into a reproducible state.
}
 \item{formula}{an object of class \code{formula} describing the model that should be fitted.
}
  \item{data}{a data frame containing the variables in the model. Each row of \code{data}
     must correspond to an observation, and each column to a binary variable (coded by 0 and 1) 
     or a factor (for details, see \code{recdom}) except for the column comprising
     the response, where no missing values are allowed in \code{data}. For a description of the specification
     of the response, see \code{y}.
}
  \item{recdom}{a logical value or vector of length \code{ncol(data)} comprising whether a SNP should
     be transformed into two binary dummy variables coding for a recessive and a dominant effect.
     If \code{recdom} is \code{TRUE} (and a logical value), then all factors/variables with three levels will be coded by two dummy
     variables as described in \code{\link[logicFS]{make.snp.dummy}}. Each level of each of the other factors 
     (also factors specifying a SNP that shows only two genotypes) is coded by one indicator variable. 
     If \code{recdom} is\code{FALSE} (and a logical value),
     each level of each factor is coded by an indicator variable. If \code{recdom} is a logical vector,
     all factors corresponding to an entry in \code{recdom} that is \code{TRUE} are assumed to be SNPs
     and transformed into two binary variables as described above. All variables corresponding
     to entries of \code{recdom} that are \code{TRUE} (no matter whether \code{recdom} is a vector or a value)
     must be coded either by the integers 1 (coding for the homozygous reference genotype), 2 (heterozygous), 
     and 3 (homozygous variant), or alternatively by the number of minor alleles, i.e. 0, 1, and 2, where
     no mixing of the two coding schemes is allowed. Thus, it is not allowed that some SNPs are coded by
     1, 2, and 3, and others are coded by 0, 1, and 2.
}
  \item{\dots}{for the \code{trioPrepare} and the \code{formula} method, optional parameters to be passed to 
    the low level function \code{trioFS.default}, i.e. all arguments of \code{trioFS.default} except for
    \code{x} and \code{y}. Otherwise, ignored.
}
}

\value{An object of class \code{trioFS} consisting of
\item{vim}{a numeric vector containing the values of the importance measure for the found interactions,}
\item{prop}{a numeric vector consisting of the percentage of models that contain the respective found interactions,}
\item{primes}{a character vector naming the found interactions,}
\item{param}{a list of parameters used in the trioFS analysis, i.e. \code{B}, \code{nleaves}, and the sampling method,}
\item{mat.imp}{if \code{addMatImp = TRUE}, a matrix containing the \code{B} improvements for each found interaction,}
\item{logreg.model}{if \code{addModel = TRUE}, the \code{B} trio logic regression models,}
\item{inbagg}{if \code{addModel = TRUE}, a list of length \code{B} in which each object specifies the trios used to
fit the corresponding trio logic regression model.}
}

\references{
Li, Q., Fallin, M.D., Louis, T.A., Lasseter, V.K., McGrath, J.A., Avramopoulos, D., Wolyniec, P.S., Valle, D., 
Liang, K.Y., Pulver, A.E., and Ruczinski, I. (2010). Detection of SNP-SNP Interactions in Trios of Parents 
with Schizophrenic Children. \emph{Genetic Epidemiology}, 34, 396-406.

Schwender, H., Bowers, K., Fallin, M.D., and Ruczinski, I. (2011). Importance Measures for Epistatic Interactions#
in Case-Parent Trios. \emph{Annals of Human Genetics}, 75, 122-132.
}
\author{Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{trioLR}}, \code{\link{print.trioFS}}, \code{\link{trio.prepare}}
}

\examples{
# Load the simulated data.
data(trio.data)

# Prepare the data in trio.ped1 for a trioFS analysis 
# by first calling
trio.tmp <- trio.check(dat = trio.ped1)

# and then applying
set.seed(123456)
trio.bin <- trio.prepare(trio.dat=trio.tmp, blocks=c(1,4,2,3))

# where we here assume the block structure to be
# c(1, 4, 2, 3), which means that the first LD "block"
# only consists of the first SNP, the second LD block
# consists of the following four SNPs in trio.bin,
# the third block of the following two SNPs,
# and the last block of the last three SNPs.
# set.seed() is specified to make the results reproducible.

# For the application of trioFS, some parameters of trio 
# logic regression are changed to make the following example faster.
my.control <- lrControl(start=1, end=-3, iter=1000, output=-4)

# Please note typically you should consider much more
# than 1000 iterations (usually, at least a few hundred
# thousand).

# TrioFS can then be applied to the trio data in trio.ped1 by
fs.out <- trioFS(trio.bin, control=my.control, rand=9876543)

# where we specify rand just to make the results reproducible.
}


\keyword{tree}
\keyword{multivariate}
\keyword{regression}
