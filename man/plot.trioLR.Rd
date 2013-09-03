\name{plot.trioLR}
\alias{plot.trioLR}
\title{Plotting for trioLR Objects}
\description{
Plots the logic trees or information on the visited models generated in a the trio logic regression analyis with \code{trioLR}. 
}

\usage{
\method{plot}{trioLR}(x, whichTree = NA, freqType = 1, useNames = FALSE, 
        addStats = TRUE, digits = 3, main = NULL, cexOper=1.5, 
		    cexLeaf=1.5, sizeLeaf=7, cexPar=1.3, ...)
}

\arguments{
  \item{x}{an object of class \code{trioLR}, i.e.\ the output of \code{\link{trioLR}}.
}
  \item{whichTree}{positive integer specifying the model for which the logic tree should be plotted
    when several trio logic regression models with different maximum numbers of leaves have been fitted.
    Ignored if just one model has been fitted using simulated annealing or MCMC has been employed to perform
    a Trio Logic Regression.} 
  \item{freqType}{positive integer between 1 and 3 specifying which statistics from the MC Trio Logic Regression
    analysis should be plotted. If \code{freqType = 1}, then for each variable, the percentage of models visited 
    (after the burn-in) in the MCMC chain that contain this variable will be plotted. If \code{freqType = 2}, then
    for each pair of variables, this percentage will be shown. If \code{freqType = 3}, then for each pair of variables,
    the observed-to-expected ratio for being jointly in the models will be plotted. Ignored if simulated annealing
    or the greedy algorithm was used in the application of \code{\link{trioLR}}.
}
 \item{useNames}{should the names of the variables be used in the plots? If \code{FALSE}, the index of the column
    is shown.}
 \item{addStats}{should the coefficient in the trio logic regression model and the score for the fitted model
    be shown in the plot? Ignored if MCMC has been used in \code{\link{trioLR}}.}
 \item{digits}{number of digits used in the presentation of the coefficient and score (see \code{addStats}). Ignored
    if \code{addStats = FALSE} or MCMC has been used in \code{\link{trioLR}}.}
 \item{main}{character string specifying the title that should be added to the plot. If \code{NULL}, a standard title
   will be added to the plot.}
 \item{cexOper}{the relative size of the AND- and OR-operators in the plotting of the logic tree. Ignored if MCMC
   has been used in \code{\link{trioLR}}.}
 \item{cexLeaf}{the relative size of the variable names shown in the logic tree. Ignored if MCMC has been used in \code{\link{trioLR}}.}
 \item{sizeLeaf}{the relative size of the boxes representing the leaves in the logic trees. Ignored if MCMC has been
   used in \code{\link{trioLR}}.}
 \item{cexPar}{the relative size of the coefficient and the score (see \code{addStats}) when plotting the logic tree.
   Ignored if \code{addStats = FALSE} or if MCMC has been used in \code{\link{trioLR}}.}
  \item{\dots}{ignored.
}
}

\references{
Kooperberg, C. and Ruczinski, I. (2005). Identifying Interacting SNPs Using Monte Carlo Logic Regression.
\emph{Genetic Epidemiology}, 28, 157-170.

Li, Q., Fallin, M.D., Louis, T.A., Lasseter, V.K., McGrath, J.A., Avramopoulos, D., Wolyniec, P.S., Valle, D., 
Liang, K.Y., Pulver, A.E., and Ruczinski, I. (2010). Detection of SNP-SNP Interactions in Trios of Parents 
with Schizophrenic Children. \emph{Genetic Epidemiology}, 34, 396-406.

Ruczinski, I., Kooperberg, C., and LeBlanc, M.L. (2003). Logic Regression. \emph{Journal of Computational and 
Graphical Statistics}, 12, 475-511.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}, based on the \code{plot} functions
implemented by Ingo Ruczinski and Charles Kooperberg in the \code{R} package \code{LogicReg}.
}


\seealso{\code{\link{trioLR}}
}

\examples{
# Load the simulated data.
data(trio.data)

# Prepare the data in trio.ped1 for a trio logic
# regression analysis by first calling
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

# For the application of trio logic regression, some
# parameters of trio logic regression are changed
# to make the following example faster.
my.control <- lrControl(start=1, end=-3, iter=1000, output=-4)

# Please note typically you should consider much more
# than 1000 iterations (usually, at least a few hundred
# thousand).

# Trio regression can then be applied to the trio data in
# trio.ped1 by
lr.out <- trioLR(trio.bin, control=my.control, rand=9876543)

# where we specify rand just to make the results reproducible.

# The logic tree representing the logic expression found in
# the trio logic regression analysis can then be plotted by
plot(lr.out)
}


\keyword{hplot}
