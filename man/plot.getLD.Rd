\name{plot.getLD}
\alias{plot.getLD}
\alias{plot.getLDlarge}

\title{
 Plotting a getLD or getLDlarge Object
}
\description{
Plots either the pairwise r^2 or D' values computed by either \code{getLD} or \code{getLDlarge}. Can also be used to plot
the categorizations used in the procedure of Gabriel et al. (2002).
}
\usage{
\method{plot}{getLD}(x, y = "rSquare", start = 1, end = NA, squared = TRUE,
   col = NULL, xlab = "", ylab = "", cexAxis = 0.8, alpha = 0.1,
   ciLD = c(0.7, 0.98), cuRecomb = 0.9, ...)
   
\method{plot}{getLDlarge}(x, y = "rSquare", start = NA, end = NA, squared = TRUE, 
   col = NULL, xlab = "", ylab = "", cexAxis = 0.8, alpha = 0.1, 
   ciLD = c(0.7,0.98), cuRecomb = 0.9, ...)
}

\arguments{
  \item{x}{the output of \code{getLD} or \code{getLDlarge}.
}
  \item{y}{either \code{"rSquare"} (default), \code{"Dprime"}, or \code{"gabriel"} specifying
    the LD values that should be plotted.
}
  \item{start}{integer or character string specifying the index or the name of the first SNP, respectively, 
    that should be plotted, where the index corresponds to the column (or row if \code{snp.in.col = FALSE})
    in the matrix used as input in \code{getLD} or \code{getLDlarge}.
}
  \item{end}{integer or character string specifying the index or the name of the last SNP, respectively,
    that should be plotted.
}
  \item{squared}{should the \eqn{r^2}{r^2} values be plotted? If \code{FALSE}, the \eqn{r} values
    are plotted. Only considered if \code{y = "rSquare"}.
}
  \item{col}{a vector specifying the colors used in plotting of the LD values. If \code{y = "rSquare"}
    or \code{y = "Dprime"}, different levels of gray will be used by default (the darker, the higher is the LD value).
    If \code{y = "gabriel"}, strong LD is by default marked by blue fields, evidence of recombination by white color,
    and others by yellow.
}
  \item{xlab}{character string naming the label of the x-axis.
}
  \item{ylab}{character string naming the label of the y-axis.
}
  \item{cexAxis}{a numeric value specifying the relative size of the SNP names displayed at the axes of the plot.
}
  \item{alpha}{numeric value between 0 and 1. Only considered if \code{y = "gabriel"}.
     For each pair of SNPs, a two-sided 
     100 * (1 - \code{alpha})\% confidence interval of D' is computed, and used
     to specify pairs of SNPs that are either in strong LD, or show historical
     evidence of recombination (see \code{ciLD} and \code{cuRecomb}). All SNP
     pairs not falling into these two categories are specified as 'Others'.  
}
  \item{ciLD}{numeric vector consisting of two values between 0 and 1. 
     Only considered if \code{y = "gabriel"}. If the 
     lower bound of the confidence interval of D' for a SNP pair is larger than
     or equal to the first value in \code{ciLD} and the upper bound is larger
     than or equal to the second value, then this pair of SNP is considered to
     be in strong LD.
}
  \item{cuRecomb}{numeric value between 0 and 1. Only considered if \code{y = "gabriel"}.
     If the upper bound of the confidence
     interval of D' for a SNP pair is smaller than \code{cuRecomb}, then this pair
     of SNP is considered to show evidence of recombination.
} 
  \item{\dots}{further arguments of \code{\link{image}}
}
}

\references{
  Gabriel, S.B. et al. (2002). The Structure of Haplotype Blocks in the Human Genome. \emph{Science}, 296, 2225-2229.
}
\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}


\seealso{
 \code{\link{getLD}}, \code{\link{plot.LDblocks}}
}

\examples{
# Load the simulated data.
data(trio.data)

# The values of Dprime and Rsquare for each pair of SNPs 
# in LDdata can be computed by
ld.out <- getLD(LDdata)

# By default, the LD measures are returned as a vector.
# If they should be returned as a matrix, then use
ld.out2 <- getLD(LDdata, asMatrix = TRUE)

# The matrix of the Rsquare values can be plotted by
plot(ld.out)

# The matrix of the Dprime values can be plotted by
plot(ld.out, "Dprime")
}

\keyword{hplot}
