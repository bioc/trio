\name{print.trioLR}
\alias{print.trioLR}
\title{Printing of trioLR Objects}
\description{
Prints information on the trio logic regression model(s) fitted with \code{\link{trioLR}}. 
}

\usage{
\method{print}{trioLR}(x, asDNF=FALSE, posBeta=FALSE, digits = 3, ...)
}

\arguments{
  \item{x}{an object of class \code{trioLR}, i.e.\ the output of \code{\link{trioLR}}.
}
  \item{asDNF}{should the disjunctive normal form of the logic expression represented by the logic tree be printed?
    If \code{FALSE}, the logic expression is printed as found by the search algorithm in trio logic regression.
    An advantage of the disjunctive normal form representation is that the interactions comprised by the logic
    expression are given by the AND-combinations in the disjunctive normal form. Note that not necessarily the
    minimum disjunctive normal form is printed so that all interactions comprised by the model are shown, even if
    some of the interactions are redundant for the evaluating the logic tree.
}
  \item{posBeta}{should the disjunctive normal form be determined as if the sign of the coefficient in trio logic
    regression model is positive? If \code{FALSE}, the sign is ignored when transforming the logic tree into
    its disjunctive normal form. If \code{TRUE} and the coefficient is negative, the complement of the logic expression is
    transformed into its disjunctive normal form and the coefficient is multiplied by -1. Ignored if \code{asDNF = FALSE}
    or the fitted logic tree only contains one leaf.
}
  \item{digits}{number of digits used in the printing of the score and the parameter estimate of the fitted trio
    logic regression model(s).
} 
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

# The output of trioLR can then be displayed by
lr.out

# This output shows the detected logic expression. If this
# expression should be displayed in disjunctive normal form,
# then this can be done by
print(lr.out, asDNF = TRUE)
}



\keyword{print}
