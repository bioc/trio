\name{poly4root}
\alias{poly4root}
\alias{poly4rootMat}

\title{
Roots of a Fourth Degree Polynomial
}
\description{
While \code{poly4root} computes the (real-valued) roots of a polynomial of fourth degree, \code{poly4rootMat} can
be applied to several polynomials of fourh degree at once by assuming that each row the input matrix contains the
coefficients for one of the polynomials.
}
\usage{
poly4root(a)

poly4rootMat(amat)
}

\arguments{
  \item{a}{a numeric vector of length five specifying the coefficients of the polynomial \code{a[1]}*\eqn{x^4}{x^4} + 
    \code{a[2]}*\eqn{x^3}{x^3} + \code{a[3]}*\eqn{x^2}{x^2} + \code{a[4]}*\eqn{x}{x} + a[5].
}
  \item{amat}{a numeric matrix with five columns in which each row contains the five coefficients of a polynomial of
    fourth degree.} 
}


\value{For \code{poly4root}, a vector containing the real-valued roots of the polynomial. For \code{poly4rootMat},
  a matrix with four columns in which each row contains the real-valued roots of the corresponding polynomial. If
  a polynomial has less than four real-valued roots, the remaining entries in the corresponding row are set to \code{NA}.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}
}

\examples{
# The roots of 
# 2 * x^4 + 3 * x^3 - x^2 + 5 * x^1 - 4 
# can be determined by
poly4root(c(2, 3, -1, 5, -4))
}  

\keyword{math}
\keyword{array}
