\name{lrControl}
\alias{lrControl}

\title{
Control Parameters for Trio Logic Regression
}
\description{
Specifies the control parameters for the search algorithms (i.e. either simulated annealing or MCMC)
and the logic tree considered when fitting a trio logic regression model. 
}
\usage{
lrControl(start = 0, end = 0, iter = 0, earlyout = 0, update = 0, 
   treesize = 8, opers = 1, minmass = 0, nburn = 1000, hyperpars = 0, 
   output = 4)
}

\arguments{
  \item{start}{a numeric value specifying the upper temperature (on log10 scale) used as start temperature in simulated annealing.
    Must be larger than \code{end}. If both \code{start = 0} and \code{end = 0}, these temperatures will be chosen automatically
    (which is not the optimal way to specify these parameters).
}
  \item{end}{a numeric value specifying the lowest temperature (on log10 scale) used in simulated annealing. Must be smaller than
    \code{start}. If both \code{start = 0} and \code{end = 0}, these temperatures will be chosen automatically (which is not
    the optimal way to specify these parameters).
}
  \item{iter}{the number of iterations used in the (stochastic) search for the best trio logic regression model, i.e. either in simulated
    annealing (if the argument \code{search} in \code{\link{trioLR}} or \code{\link{trioFS}} is set to \code{"sa"}) or
    in MCMC (if \code{search = "mcmc"}). If \code{iter = 0}, \code{iter} will be chosen automatically (similar to \code{start} and \code{end})
    when simulated annealing is used, and will be set to \code{iter = 50000} when MCMC is employed.
}
  \item{earlyout}{a non-negative integer providing an option to end the search before all \code{iter} iterations in simulated annealing are considered.
    If during five consecutive blocks of \code{earlyout} iterations, 10 or fewer moves proposed in simulated annealing are accepted in
    each of the blocks, then the search will terminate. Can help to stop the search earlier, when there is no progress in the search
    anymore. By default, all \code{iter} iterations are considered. 
}
  \item{update}{the number of iterations in simulated annealing or MCMC after which statistics for the current trio logic regression model are
    displayed. This argument allows to evaluate the progress in the search for the best trio logic regression model. By default,
    no updates are shown.
}
  \item{treesize}{a positive integer specifying the maximum number of leaves allowed in the logic tree of a trio logic regression model.
}
  \item{opers}{either 1, 2, or 3 specifying if both the AND and the OR operator (\code{opers = 1}), or only the AND operator
    (\code{opers = 2}), or only the OR operator (\code{opers = 3}) is considered when building the logic tree.
}
  \item{minmass}{a non-negative integer specifying the number of cases and pseudo-controls for which the logic expression (i.e. the logic tree)
    needs to be 1 or for which the logic expression needs to be 0 to be considered as a logic tree in the trio logic regression model.
    By default, \code{minmass} is either set to 20\% of the trios or to 15, whatever is less.
}
  \item{nburn}{number of initial iterations in MCMC considered as burn-in MC trio logic regression, and therefore, ignored when computing
    the summaries.
}
  \item{hyperpars}{a numeric value specifying the hyperparameter for the prior on the model size when performing a MC trio logic regression.
    More exactly, \code{hyperpars} is assumed to be \eqn{log(P(size = k) / P(size = k+1))}, where \eqn{P} is the prior on the model size.
}
  \item{output}{a value specifying which statistics are returned in an MCMC trio logic regression analysis. If \code{output > 0}, then all
    fitted models are saved in a text file called "triolrlisting.tmp" in the current working directory. By setting \code{output < 0}, this
    can be avoided. If \code{abs(output) > 1}, bivariate statistics are gathered. If \code{abs(output) > 2}, trivariate statistics are
    gathered. Otherwise, only univariate statistics are determined. 
}
}
\details{More details on the different control parameters and their specification can be found on the help pages of the functions
  \code{logreg.anneal.control}, \code{logreg.tree.control}, and \code{logreg.mc.control} for the different
  types of control parameters available in the \code{R} package \code{LogicReg} for a standard logic regressions.  
}
\value{A list containing all required control parameters.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\examples{
# The default values for the parameters in trio logic regression
# can be specified by
myControl <- lrControl()

# If the starting temperature of Simulated Annealing should be set
# to 100 and the lowest temperature to 0.001, then this can be done by
myControl2 <- lrControl(start = 2, end = -3)
}

\keyword{utilities}
